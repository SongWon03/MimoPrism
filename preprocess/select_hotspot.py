#!/usr/bin/env python
"""
<select_hotspot.py>

Input : processed_<PDB>.pdb (output of preprocess.py)
1) base      : no hotspot (control)
2) cbeta     : antibody residues whose Cβ is within --cbeta-cutoff Å of an antigen Cβ, nearest first
3) contact   : antibody residues with a heavy atom within --contact-cutoff Å of the antigen,
               ranked by ΔSASA (area buried by the antigen) or by distance
4) parasurf  : extract H/L → run ParaSurf → residues ranked by paratope probability
5) write design_names / hotspots (+ PEPTIDE_CHAIN_ID, CONTIGS_RECEPTOR) into the config JSON

Example:
  python select_hotspot.py 6AXK_Fv --ab H L --ag P --ks 6
  python select_hotspot.py 6AXK_Fv --ab H L --ag P --ks 6 --sources contact parasurf
  python select_hotspot.py 6AXK_Fv --ab H L --ag P --dry-run        # print lists only, config untouched
"""
import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# Layout: <MimoPrism>/preprocess/select_hotspot.py
SCRIPT_DIR = Path(__file__).resolve().parent               # preprocess/
PROJECT_DIR = SCRIPT_DIR.parent                            # MimoPrism/  (input/, tool/, *.json)
DEFAULT_DIR = PROJECT_DIR / "input"                        # processed_<PDB>.pdb
DEFAULT_TOOL_DIR = PROJECT_DIR / "tool"                    # external tools (ParaSurf, ...)
DEFAULT_ENV_JSON = PROJECT_DIR / "mimoprism_env.json"
DEFAULT_CONFIG = PROJECT_DIR / "mimoprism_config.json"
SOURCES = ("base", "cbeta", "contact", "parasurf")


# ── environment (mimoprism_env.json) ────────────────────────────────────────
def load_env(path):
    """Flatten mimoprism_env.json into {KEY: value}. Exported environment variables win,
    then the JSON, then the caller's fallback. Same priority as mp_load_config.py --soft."""
    cfg = {}
    if path and Path(path).is_file():
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        for section, params in raw.items():
            if section.startswith("_") or not isinstance(params, dict):
                continue
            for k, v in params.items():
                if k.startswith("_") or v is None:
                    continue
                cfg[k] = os.path.expanduser(os.path.expandvars(v)) if isinstance(v, str) else v
    for k in list(cfg):
        if k in os.environ:
            cfg[k] = os.environ[k]
    return cfg


def env_get(env, key, fallback):
    return os.environ.get(key) or env.get(key) or fallback


# ── shared helpers ──────────────────────────────────────────────────────────
def reskey(res):
    """(chain, resnum, icode). Insertion codes matter: CDR-L1 27A-27D would collapse without them."""
    return (res.get_parent().id, res.id[1], res.id[2])


def label(res):
    """RFdiffusion notation: chain + number + insertion code (H31, L27A)."""
    ch, num, ic = reskey(res)
    return f"{ch}{num}{ic.strip()}"


def load_model(path, drop_het=True):
    """PDB → first model. drop_het removes water/ligands (SASA depends on what is left)."""
    from Bio.PDB import PDBParser
    st = PDBParser(QUIET=True).get_structure("x", str(path))
    model = st[0]
    if drop_het:
        for ch in model:
            for r in [r for r in ch if r.id[0] != " "]:
                ch.detach_child(r.id)
    return model


def cb_or_ca(res):
    return res["CB"] if "CB" in res else (res["CA"] if "CA" in res else None)


def fmt_list(labels):
    """['H31','L92',...]  (no spaces — the format mimoprism_config.json / RFdiffusion expect)."""
    return "[" + ",".join(f"'{x}'" for x in labels) + "]"




def select(scored, coords, top, min_sep):
    """Top `top` entries of a score-sorted list. With min_sep > 0 a residue closer than min_sep Å
    (Cβ) to an already chosen one is skipped, spreading picks over the paratope."""
    picked = []
    for key, lab, val in scored:
        if top and len(picked) >= top:
            break
        if min_sep > 0 and key in coords:
            p = coords[key]
            if any(q in coords and math.dist(p, coords[q]) < min_sep for q, _, _ in picked):
                continue
        picked.append((key, lab, val))
    return picked


def report(name, picked, n_total, value_name, value_fmt, unit=""):
    print(f"  [{name}] candidates {n_total} → selected {len(picked)}  (ranked by {value_name})")
    for _k, lab, v in picked:
        print(f"    {lab:<8} {v:{value_fmt}} {unit}".rstrip())
    print(f"  [{name}] {fmt_list([lab for _k, lab, _v in picked])}")


# ── native complex (cbeta / contact) ────────────────────────────────────────
def split_ab_ag(pdb, ab, ag):
    """Complex → (model, antibody residues, antigen residues, Cβ coordinate table)."""
    model = load_model(pdb)
    ab, ag = set(ab), set(ag)
    ab_res = [r for ch in model if ch.id in ab for r in ch]
    ag_res = [r for ch in model if ch.id in ag for r in ch]
    present = sorted(ch.id for ch in model)
    hint = ("\n        chain IDs do not match: either pass --ab/--ag, or re-run preprocess.py — it renames "
            "chains to H/L (antibody) and P (antigen), so an older processed PDB still carries the "
            "original IDs")
    if not ab_res:
        raise SystemExit(f"[FATAL] no residues in antibody chains {sorted(ab)} (chains in PDB: {present}){hint}")
    if not ag_res:
        raise SystemExit(f"[FATAL] no residues in antigen chains {sorted(ag)} (chains in PDB: {present}){hint}")
    coords = {reskey(r): tuple(map(float, a.coord)) for r in ab_res
              if (a := cb_or_ca(r)) is not None}
    return model, ab_res, ag_res, coords




def run_cbeta(pdb, ab, ag, cutoff):
    """RFdiffusion definition: antibody Cβ (Cα for Gly) within cutoff of any antigen Cβ, nearest first.
    No --top / --min-sep: this reproduces the paper's method 1 as is."""
    _model, ab_res, ag_res, _coords = split_ab_ag(pdb, ab, ag)
    ag_atoms = [a for r in ag_res if (a := cb_or_ca(r)) is not None]
    hits = []
    for r in ab_res:
        a = cb_or_ca(r)
        if a is None:
            continue
        dmin = min(a - g for g in ag_atoms)
        if dmin <= cutoff:
            hits.append((reskey(r), label(r), dmin))
    hits.sort(key=lambda x: x[2])
    report("cbeta", hits, len(hits), f"Cβ distance ≤ {cutoff:g} Å, ascending", "5.2f", "Å")
    return [lab for _k, lab, _v in hits]


def dsasa_table(pdb, ab, ag):
    """{reskey: (ΔSASA, SASA_free)} — Shrake-Rupley (probe 1.4 Å) on the complex, then on the same
    coordinates with the antigen chains removed. ΔSASA = SASA_free − SASA_complex."""
    from Bio.PDB.SASA import ShrakeRupley
    sr = ShrakeRupley()
    m_c = load_model(pdb)
    sr.compute(m_c, level="R")
    s_c = {reskey(r): r.sasa for ch in m_c if ch.id in ab for r in ch}
    m_f = load_model(pdb)
    for cid in [ch.id for ch in m_f if ch.id in ag]:
        m_f.detach_child(cid)
    sr.compute(m_f, level="R")
    s_f = {reskey(r): r.sasa for ch in m_f if ch.id in ab for r in ch}
    return {k: (v - s_c.get(k, 0.0), v) for k, v in s_f.items()}


def run_contact(pdb, ab, ag, cutoff, rank, top, min_sep):
    """Antibody residues with any heavy atom within cutoff of the antigen, ranked by ΔSASA or distance."""
    _model, ab_res, ag_res, coords = split_ab_ag(pdb, ab, ag)
    ag_atoms = [a for r in ag_res for a in r]
    hits = []
    for r in ab_res:
        dmin = min((a - g for a in r for g in ag_atoms), default=1e9)
        if dmin <= cutoff:
            hits.append((reskey(r), label(r), dmin))
    if not hits:
        raise SystemExit(f"[FATAL] no antibody residue within {cutoff} Å (heavy atom) of the antigen")
    if rank == "dist":
        hits.sort(key=lambda x: x[2])
        scored, vname, vfmt, unit = hits, "min heavy-atom distance, ascending", "5.2f", "Å"
    else:
        d = dsasa_table(pdb, set(ab), set(ag))
        scored = [(k, l, d[k][0]) for k, l, _v in hits if k in d and d[k][0] > 0]
        scored.sort(key=lambda x: -x[2])
        vname, vfmt, unit = "ΔSASA, descending", "6.1f", "Å²"
        print(f"  [contact] residues with ΔSASA > 0: {len(scored)}/{len(hits)}, "
              f"antibody-side ΔSASA total {sum(v for _k, _l, v in scored):.0f} Å²")
    picked = select(scored, coords, top, min_sep)
    report("contact", picked, len(scored), vname, vfmt, unit)
    return [lab for _k, lab, _v in picked]


# ── ParaSurf ────────────────────────────────────────────────────────────────
def extract_chains(pdb, out, keep):
    """Write only the given chains (standard residues) — ParaSurf input."""
    from Bio.PDB import PDBParser, PDBIO, Select
    keep = set(keep)

    class Sel(Select):
        def accept_chain(self, chain):
            return 1 if chain.id in keep else 0

        def accept_residue(self, residue):
            return 1 if residue.id[0] == " " else 0

    st = PDBParser(QUIET=True).get_structure("c", str(pdb))
    io = PDBIO(); io.set_structure(st)
    io.save(str(out), Sel())
    print(f"  [parasurf] antibody-only PDB → {out}  (chains {sorted(keep)})")


def run_parasurf_predict(hl_pdb, ps):
    """Run ParaSurf blind_predict.py in its own interpreter. Returns (pred_pdb, pocket_pdb)."""
    script = Path(ps["dir"]) / ps["script"]
    if not script.is_file():
        raise SystemExit(f"[FATAL] ParaSurf script not found: {script}\n"
                         "        clone ParaSurf into <MimoPrism>/tool/ParaSurf (or pass --parasurf-dir), "
                         "or run without it: --sources base cbeta contact")
    env = os.environ.copy()
    if ps["dms_bin_dir"]:
        env["PATH"] = f"{ps['dms_bin_dir']}{os.pathsep}{env.get('PATH', '')}"
    if not shutil.which("dms", path=env["PATH"]):
        raise SystemExit("[FATAL] 'dms' (surface generation) not found on PATH — ParaSurf cannot build the surface.\n"
                         "        Find it with `which dms` in the working ParaSurf env and pass --dms-bin-dir.")
    # dms looks for the atom radius table in the cwd first, then in <LIBDIR>/dms/radii.
    if not any(Path(p).is_file() for p in (Path(ps["dir"]) / "radii", Path("/usr/local/lib/dms/radii"),
                                           Path.home() / ".local/lib/dms/radii")):
        print("  [warn] dms radius table (radii) not found in the usual places; if dms fails, copy "
              f"{ps['dir']}/dms/radii.proto to {ps['dir']}/radii", file=sys.stderr)
    # ParaSurf's README exports PYTHONPATH in ~/.bashrc; non-interactive shells never read it.
    env["PYTHONPATH"] = ps["dir"] + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    cmd = [ps["python"], ps["script"], "--receptor", str(hl_pdb.resolve()), "--model_weights", ps["weights"]]
    print(f"  [parasurf] running: {' '.join(cmd)}  (cwd={ps['dir']})")
    r = subprocess.run(cmd, cwd=ps["dir"], env=env)
    if r.returncode != 0:
        raise SystemExit(f"[FATAL] ParaSurf exited with code {r.returncode}")

    hl_dir, stem = hl_pdb.parent, hl_pdb.stem
    pred = next((c for c in (hl_dir / stem / "HL_pred.pdb", hl_dir / stem / f"{stem}_pred.pdb",
                             hl_dir / f"{stem}_pred.pdb") if c.is_file()), None)
    if pred is None:
        raise SystemExit("[FATAL] ParaSurf prediction PDB not found — check where ParaSurf wrote its output "
                         f"(looked under {hl_dir / stem} and {hl_dir})")
    pockets = sorted((hl_dir / stem).glob("pocket*.pdb")) + sorted(hl_dir.glob(f"{stem}*pocket*.pdb"))
    pocket = pockets[0] if pockets else None
    print(f"  [parasurf] pred={pred}  pocket={pocket or 'none (ranking over all residues)'}")
    return pred, pocket


def pdb_reskey(line):
    """PDB ATOM line → (chain, resnum, icode)."""
    return (line[21], int(line[22:26]), line[26])


def run_parasurf(pred, pocket, top, min_sep):
    """Residue probability = max B-factor of its atoms in the prediction PDB, descending.
    Probability only — no area weighting; --min-sep keeps picks from piling into one CDR loop."""
    pocket_keys = None
    if pocket:
        pocket_keys = {pdb_reskey(l) for l in open(pocket) if l.startswith(("ATOM", "HETATM"))}
    prob, coords = {}, {}
    for l in open(pred):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        k = pdb_reskey(l)
        try:
            prob[k] = max(prob.get(k, -1.0), float(l[60:66]))
        except ValueError:
            continue
        atom = l[12:16].strip()
        if atom == "CB" or (atom == "CA" and k not in coords):
            coords[k] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    scored = [(k, f"{k[0]}{k[1]}{k[2].strip()}", s) for k, s in prob.items()
              if pocket_keys is None or k in pocket_keys]
    scored.sort(key=lambda x: -x[2])
    picked = select(scored, coords, top, min_sep)
    report("parasurf", picked, len(scored), "ParaSurf probability, descending", ".3f")
    return [lab for _k, lab, _v in picked]


# ── config ──────────────────────────────────────────────────────────────────
def build_designs(pdb, lists, ks, suffix, mode):
    """design_names / hotspots. lists = [(source, labels)] in order; base has an empty list.
    An empty non-base source is skipped; K larger than the list is skipped too (a silently
    truncated Hot7 that is really Hot5 would corrupt the analysis)."""
    names, hots = [], []
    for source, lst in lists:
        if source == "base":
            names.append(f"{pdb}_{mode}_base_Hot0_{suffix}")
            hots.append("[]")
            continue
        if not lst:
            print(f"  [config] {source}: empty list → no design")
            continue
        for k in ks:
            if k > len(lst):
                print(f"  [config] {source}: K={k} > {len(lst)} candidates → skipped")
                continue
            names.append(f"{pdb}_{mode}_{source}_Hot{k}_{suffix}")
            hots.append(fmt_list(lst[:k]))
    return names, hots


def chain_lengths(pdb, chains):
    """Standard residues per chain in the processed PDB → CONTIGS_RECEPTOR."""
    model = load_model(pdb)
    return {ch.id: sum(1 for r in ch if r.id[0] == " ") for ch in model if ch.id in set(chains)}


def load_config(path):
    """Config JSON → dict. A missing or empty file is a fresh config ({}); anything else that
    is not valid JSON is an error (never overwrite a hand-edited file we cannot parse)."""
    path = Path(path)
    if not path.is_file() or not path.read_text(encoding="utf-8-sig").strip():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8-sig"))
    except json.JSONDecodeError as e:
        raise SystemExit(f"[FATAL] config is not valid JSON: {path}  ({e})\n"
                         "        fix the file or pass a different --config")


def write_config(path, pdb_name, names, hots, pep_chain, contigs):
    cfg = load_config(path)
    g = cfg.setdefault("Global", {})
    g["PDB_NAME"] = pdb_name
    g["design_names"] = names
    g["hotspots"] = hots
    if pep_chain:
        g["PEPTIDE_CHAIN_ID"] = pep_chain
    if contigs:
        cfg.setdefault("RFpeptide", {})["CONTIGS_RECEPTOR"] = contigs
    path.write_text(json.dumps(cfg, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def default_suffix(config_path):
    """1C{N}C1 from the config's SS_LOOP_N when it is a disulfide-cyclic campaign, else 1C8C1."""
    rf = load_config(config_path).get("RFpeptide", {})
    if rf.get("DISULFIDE_CYCLIC", True) and rf.get("SS_LOOP_N"):
        return f"1C{int(rf['SS_LOOP_N'])}C1"
    print(f"  [warn] SS_LOOP_N not found in {config_path} → suffix falls back to 1C8C1 (override with --suffix)",
          file=sys.stderr)
    return "1C8C1"


# ── main ────────────────────────────────────────────────────────────────────
def main():
    env = load_env(os.environ.get("MP_ENV", DEFAULT_ENV_JSON))
    ap = argparse.ArgumentParser(
        description="hotspot candidates (base / cbeta / contact / parasurf) → mimoprism config")
    ap.add_argument("pdb_name", help="processed name, e.g. 6AXK_Fv → <dir>/processed_6AXK_Fv.pdb")
    ap.add_argument("--dir", default=None,
                    help="folder with processed_<PDB>.pdb (default: <MimoPrism>/input — "
                         "the same folder preprocess.py writes to)")
    ap.add_argument("--ab", nargs=2, metavar=("H", "L"), default=["H", "L"], help="antibody chain IDs")
    ap.add_argument("--ag", nargs="+", default=["P"], help="antigen chain IDs")
    ap.add_argument("--pep", default=None,
                    help="PEPTIDE_CHAIN_ID written to the config (default: first --ag chain)")
    ap.add_argument("--sources", nargs="+", choices=SOURCES, default=list(SOURCES),
                    help="which hotspot sources to run (default: all four)")
    # selection
    ap.add_argument("--ks", nargs="+", type=int, default=[3, 4, 5, 6, 7], help="hotspots per design")
    ap.add_argument("--min-sep", dest="min_sep", type=float, default=0.0,
                    help="min Cβ distance (Å) between chosen hotspots, contact/parasurf only (0 = off)")
    ap.add_argument("--cbeta-cutoff", dest="cbeta_cutoff", type=float, default=10.0)
    ap.add_argument("--contact-cutoff", dest="contact_cutoff", type=float, default=4.5)
    ap.add_argument("--contact-rank", dest="contact_rank", choices=["dsasa", "dist"], default="dsasa")
    # config
    ap.add_argument("--config", default=os.environ.get("MP_CONFIG", str(DEFAULT_CONFIG)),
                    help="campaign config to update (default: <MimoPrism>/mimoprism_config.json)")
    ap.add_argument("--suffix", default=None, help="design-name suffix (default: 1C{SS_LOOP_N}C1 from the config)")
    ap.add_argument("--mode", default="cyclic", help="design-name mode token")
    ap.add_argument("--dry-run", dest="dry_run", action="store_true", help="print lists only; config untouched")
    ap.add_argument("--parasurf-dir", default=str(DEFAULT_TOOL_DIR / "ParaSurf"),
                    help="ParaSurf folder (default: <MimoPrism>/tool/ParaSurf)")
    ap.add_argument("--parasurf-script", default="blind_predict.py")
    ap.add_argument("--parasurf-weights", default="ParaSurf_best.pth")
    ap.add_argument("--parasurf-python", default=env_get(env, "PYTHON_PARASURF", env_get(env, "PYTHON", "python")))
    ap.add_argument("--dms-bin-dir", default=str(DEFAULT_TOOL_DIR / "ParaSurf" / "dms"),
                    help="folder holding the dms binary (default: <MimoPrism>/tool/ParaSurf/dms)")
    args = ap.parse_args()

    d = Path(args.dir or DEFAULT_DIR).expanduser()
    pdb = (d / f"processed_{args.pdb_name}.pdb").resolve()
    if not pdb.is_file():
        raise SystemExit(f"[FATAL] processed PDB not found: {pdb}  (run preprocess.py first, or pass --dir)")
    config_path = Path(args.config).expanduser()
    if not config_path.is_file() and not args.dry_run:
        raise SystemExit(f"[FATAL] config not found: {config_path}\n"
                         "        pass --config <mimoprism_*.json> (or export MP_CONFIG) so the existing campaign "
                         "config is updated instead of a new file being created")
    if config_path.is_file() and not load_config(config_path):
        print(f"  [warn] {config_path} is empty → a new config with Global/RFpeptide only will be written; "
              "the driver still needs the other sections", file=sys.stderr)
    suffix = args.suffix or default_suffix(config_path)
    pep_chain = args.pep or args.ag[0]
    print(f"[select_hotspot] pdb={pdb}  config={config_path}")
    print(f"[select_hotspot] ab={args.ab} ag={args.ag}  sources={args.sources}  ks={args.ks}  suffix={suffix}")

    lists = []
    for source in args.sources:
        if source == "base":
            lists.append(("base", []))
        elif source == "cbeta":
            print(f"[cbeta] Cβ ≤ {args.cbeta_cutoff:g} Å")
            lists.append(("cbeta", run_cbeta(pdb, args.ab, args.ag, args.cbeta_cutoff)))
        elif source == "contact":
            print(f"[contact] heavy atom ≤ {args.contact_cutoff:g} Å, rank={args.contact_rank}, "
                  f"min_sep={args.min_sep:g} Å")
            lists.append(("contact", run_contact(pdb, args.ab, args.ag, args.contact_cutoff,
                                                 args.contact_rank, 0, args.min_sep)))
        elif source == "parasurf":
            print(f"[parasurf] min_sep={args.min_sep:g} Å")
            hl_pdb = d / f"HL_only_{args.pdb_name}.pdb"
            extract_chains(pdb, hl_pdb, args.ab)
            ps = {"dir": str(Path(args.parasurf_dir).expanduser()), "script": args.parasurf_script,
                  "weights": args.parasurf_weights, "python": str(Path(args.parasurf_python).expanduser()),
                  "dms_bin_dir": args.dms_bin_dir}
            pred, pocket = run_parasurf_predict(hl_pdb, ps)
            lists.append(("parasurf", run_parasurf(pred, pocket, 0, args.min_sep)))

    print("[config]")
    names, hots = build_designs(args.pdb_name, lists, args.ks, suffix, args.mode)
    if not names:
        raise SystemExit("[FATAL] no design generated — config left untouched "
                         "(all sources empty, or every --ks larger than the candidate lists)")
    lengths = chain_lengths(pdb, args.ab)
    contigs = " ".join(f"{c}1-{lengths[c]}/0" for c in args.ab if lengths.get(c))
    for n, h in zip(names, hots):
        print(f"  {n:<44} {h}")
    print(f"  PEPTIDE_CHAIN_ID = {pep_chain}")
    print(f"  CONTIGS_RECEPTOR = \"{contigs}\"")
    if args.dry_run:
        print("[config] --dry-run: nothing written")
        return
    write_config(config_path, args.pdb_name, names, hots, pep_chain, contigs)
    print(f"[config] updated {config_path}  (Global.PDB_NAME / design_names / hotspots / PEPTIDE_CHAIN_ID, "
          f"RFpeptide.CONTIGS_RECEPTOR)")


if __name__ == "__main__":
    main()