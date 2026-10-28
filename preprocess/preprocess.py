#!/usr/bin/env python
"""
<preprocess.py>

1) Convert noncanonical residues into canonical residues
2) Convert fab to fv chain using ANARCI (DEFAULT: imgt)
3) Remove water and ligands
4) Renumber all residues

Example:
  python preprocess.py 6AXK --ab A B --ag E
"""
import argparse
import csv
import os
import re
import sys
from pathlib import Path

# ── 1) Noncanonical residue → parent canonical residue ──────────────────────
MOD_PARENT = {
    "MSE": "MET",   # selenomethionine
    "PCA": "GLN",   # pyroglutamate
    "SEP": "SER",   # phosphoserine
    "TPO": "THR",   # phosphothreonine
    "PTR": "TYR",   # phosphotyrosine
    "CSO": "CYS",   # S-hydroxycysteine
    "HYP": "PRO",   # 4-hydroxyproline
    "CME": "CYS",   # S,S-(2-hydroxyethylthio)cysteine
    "KCX": "LYS",   # N-carboxylysine
    "LLP": "LYS",   # lysine-pyridoxal-5'-phosphate
    "MLY": "LYS",   # N-dimethyllysine
    "SAC": "SER",   # N-acetylserine
}

# Only atoms whose name differs from the parent residue (others already match).
ATOM_RENAME = {
    "MSE": {"SE": "SD"},
    "PCA": {"OE": "OE1"},
}
BACKBONE = {"N", "CA", "C", "O", "OXT"}
# Side-chain heavy atoms of the canonical residues (backbone excluded).
SIDECHAIN = {
    "ALA": {"CB"},
    "ARG": {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
    "ASN": {"CB", "CG", "OD1", "ND2"},
    "ASP": {"CB", "CG", "OD1", "OD2"},
    "CYS": {"CB", "SG"},
    "GLN": {"CB", "CG", "CD", "OE1", "NE2"},
    "GLU": {"CB", "CG", "CD", "OE1", "OE2"},
    "GLY": set(),
    "HIS": {"CB", "CG", "ND1", "CD2", "CE1", "NE2"},
    "ILE": {"CB", "CG1", "CG2", "CD1"},
    "LEU": {"CB", "CG", "CD1", "CD2"},
    "LYS": {"CB", "CG", "CD", "CE", "NZ"},
    "MET": {"CB", "CG", "SD", "CE"},
    "PHE": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "PRO": {"CB", "CG", "CD"},
    "SER": {"CB", "OG"},
    "THR": {"CB", "OG1", "CG2"},
    "TRP": {"CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "TYR": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
    "VAL": {"CB", "CG1", "CG2"},
}
AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
WATER = {"HOH", "WAT", "H2O", "DOD"}
ANARCI_TYPE = {"H": "heavy", "K": "light", "L": "light"}
# Layout: <MimoPrism>/preprocess/preprocess.py → input folder is <MimoPrism>/input
PROJECT_DIR = Path(__file__).resolve().parent.parent
DEFAULT_DIR = PROJECT_DIR / "input"


def fmt_atom_name(name):
    """PDB columns 13-16: a 1-letter element with a name of <=3 chars starts at column 14."""
    return name if len(name) >= 4 else f" {name:<3}"


def guess_element(name):
    c = name[0]
    return c if c in "NCOSP" else " "


def convert_nonstd(lines):
    """Replace noncanonical residues. Returns (lines, report dict, dropped atom serials)."""
    out, dropped_serials = [], set()
    # "A1(PCA→GLN)" → {"hydrogens": n, "extra": [...heavy atoms not in parent...],
    #                  "missing": [...parent atoms the residue never had...]}
    report = {}
    seen_atoms = {}
    for l in lines:
        if not l.startswith(("ATOM", "HETATM")) or len(l) < 27:
            out.append(l)
            continue
        rn = l[17:20].strip().upper()
        if rn not in MOD_PARENT:
            out.append(l)
            continue
        parent = MOD_PARENT[rn]
        key = f"{l[21]}{l[22:27].strip()}({rn}→{parent})"
        rep = report.setdefault(key, {"hydrogens": 0, "extra": [], "missing": []})
        seen_atoms.setdefault(key, set())
        name = l[12:16].strip()
        name = ATOM_RENAME.get(rn, {}).get(name, name)
        if name not in BACKBONE and name not in SIDECHAIN[parent]:   # drop hydrogens / extra atoms
            if name.lstrip("0123456789")[:1] in ("H", "D"):   # hydrogen / deuterium
                rep["hydrogens"] += 1
            else:
                rep["extra"].append(l[12:16].strip())
            try:
                dropped_serials.add(int(l[6:11]))
            except ValueError:
                pass
            continue
        seen_atoms[key].add(name)
        new = ("ATOM  " + l[6:12] + fmt_atom_name(name) + l[16]
               + f"{parent:>3}" + l[20:])
        if len(new) >= 78:                      # fix element column (e.g. SE → S)
            new = new[:76] + f"{guess_element(name):>2}" + new[78:]
        out.append(new)
    for key, names in seen_atoms.items():
        parent = key.split("→")[1].rstrip(")")
        report[key]["missing"] = sorted((SIDECHAIN[parent] | {"N", "CA", "C", "O"}) - names)
    return out, report, dropped_serials


def clean_records(lines, dropped_serials, converted_names):
    """Drop CONECT records pointing at deleted atoms and HET/HETNAM/FORMUL/MODRES
    records of converted residues."""
    out = []
    for l in lines:
        if l.startswith("CONECT"):
            fields = [l[i:i + 5].strip() for i in range(6, min(len(l.rstrip()), 31), 5)]
            serials = [int(f) for f in fields if f.isdigit()]
            if any(s in dropped_serials for s in serials):
                continue
        elif l.startswith(("HET   ", "HETNAM", "HETSYN", "FORMUL", "MODRES")):
            if any(n in l[:40] for n in converted_names):
                continue
        out.append(l)
    return out


def unknown_hetatm(lines):
    """HETATM residue names not in the table (ligands, or amino-acid derivatives not listed yet)."""
    return sorted({
        l[17:20].strip().upper() for l in lines
        if l.startswith("HETATM") and len(l) >= 20
        and l[17:20].strip().upper() not in set(SIDECHAIN) | set(MOD_PARENT) | WATER
    })


# ── 2) Fab → Fv ─────────────────────────────────────────────────────────────
def anarci_domain(seq, scheme, label):
    """Locate the V-domain boundary. Returns (cut, chain_type); cut = residues to keep."""
    try:
        from anarci import run_anarci
    except Exception as e:
        raise SystemExit(
            "[FATAL] cannot import ANARCI.\n"
            "  install: conda install -c bioconda anarci   (or: pip install anarci)\n"
            f"  cause: {e}")
    try:
        out = run_anarci([("q", seq)], scheme=scheme)
    except Exception as e:
        raise SystemExit(f"[FATAL] ANARCI failed ({label}): {e}")
    if len(out) == 4:                          # ANARCI returns a 3- or 4-tuple depending on version
        _, numbering, details, _ = out
    elif len(out) == 3:
        numbering, details, _ = out
    else:
        raise SystemExit(f"[FATAL] unrecognized ANARCI return format ({label}): len={len(out)}")
    if not numbering or not numbering[0]:
        raise SystemExit(
            f"[FATAL] ANARCI found no V-domain in {label}.\n"
            f"  check that the chain is an antibody variable domain and the chain ID is right.\n"
            f"  sequence start: {seq[:60]}")
    dom = numbering[0][0]                      # (numbering, start_idx, end_idx)
    end = int(dom[2])
    ctype = None
    try:
        ctype = (details[0][0] or {}).get("chain_type")
    except Exception:
        pass
    return end + 1, ctype


def read_chains(lines):
    """chain → {"keys": [resSeq+iCode ...], "res": {key: resName}, "idx": {key: [line idx]},
    "std": {key: is ATOM record}}. Water is skipped. Step 1 already ran, so MSE etc.
    are ATOM records with canonical names here."""
    chains = {}
    for i, l in enumerate(lines):
        if not l.startswith(("ATOM", "HETATM")) or len(l) < 27:
            continue
        rn = l[17:20].strip().upper()
        if rn in WATER:
            continue
        ch, key = l[21], l[22:27]
        c = chains.setdefault(ch, {"keys": [], "res": {}, "idx": {}, "std": {}})
        if key not in c["res"]:
            c["keys"].append(key)
            c["res"][key] = rn
            c["idx"][key] = []
            c["std"][key] = l.startswith("ATOM")
        c["idx"][key].append(i)
    return chains


def std_keys(c):
    """Residues fed to ANARCI = ATOM records only; HETATM attached to the chain (e.g. glycans) excluded."""
    return [k for k in c["keys"] if c["std"][k]]


def chain_seq(c):
    return "".join(AA3TO1.get(c["res"][k], "X") for k in std_keys(c))


def assign_hl(ab, info):
    c1, c2 = ab
    t1, t2 = ANARCI_TYPE.get(info[c1]["type"]), ANARCI_TYPE.get(info[c2]["type"])
    if t1 == "heavy" and t2 == "light":
        return c1, c2
    if t1 == "light" and t2 == "heavy":
        print(f"  [info] heavy/light swapped per ANARCI chain type: H={c2}, L={c1}")
        return c2, c1
    print(f"  [warn] ambiguous chain_type ({info[c1]['type']}, {info[c2]['type']}) — "
          f"using input order H={c1}, L={c2}", file=sys.stderr)
    return c1, c2


def cut_fv(lines, ab, ag, scheme, ab_rename, ag_rename):
    """Keep antibody chains up to the V-domain, antigen chains in full, drop every other chain.
    Chain IDs are rewritten: heavy/light → ab_rename=(H_id, L_id), antigens → ag_rename (same order as ag).
    Returns (lines, original H id, original L id, rename_map {old: new})."""
    chains = read_chains(lines)
    for ch in list(ab) + list(ag):
        if ch not in chains:
            raise SystemExit(f"[FATAL] chain '{ch}' not in PDB (available: {sorted(chains)})")
    info = {}
    for ch in ab:
        seq = chain_seq(chains[ch])
        cut, ctype = anarci_domain(seq, scheme, f"PDB chain {ch}")
        info[ch] = {"cut": cut, "type": ctype, "len": len(seq)}
        keys = std_keys(chains[ch])
        print(f"  [PDB] chain {ch}: {len(seq)} → {cut} residues "
              f"(res {keys[0][:4].strip()}..{keys[cut - 1][:4].strip()})  chain_type={ctype}")
    H, L = assign_hl(ab, info)

    keep_idx = set()
    for ch in [H, L] + list(ag):
        c = chains[ch]
        if ch in info:
            keep = set(std_keys(c)[:info[ch]["cut"]])
            # HETATM attached to the chain (ligand/glycan) is kept here and handled by step 3 flags.
            keep |= {k for k in c["keys"] if not c["std"][k]}
        else:
            keep = set(c["keys"])
        for k in c["keys"]:
            if k in keep:
                keep_idx.update(c["idx"][k])

    rename_map = {H: ab_rename[0], L: ab_rename[1], **dict(zip(ag, ag_rename))}
    written, prev = [], None
    for i, l in enumerate(lines):
        if not l.startswith(("ATOM", "HETATM")):
            if l.startswith(("TER", "END", "MODEL", "ENDMDL")):
                continue
            written.append(l)                    # keep header / REMARK records
            continue
        if i not in keep_idx:
            continue
        newch = rename_map.get(l[21], l[21])
        if prev is not None and newch != prev:
            written.append("TER\n")
        prev = newch
        written.append(l[:21] + newch + l[22:])
    written.append("TER\nEND\n")
    return written, H, L, rename_map


# ── FASTA ───────────────────────────────────────────────────────────────────
def parse_fasta(path):
    header, buf, out = None, [], []
    for line in Path(path).expanduser().read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                out.append((header, "".join(buf)))
            header, buf = line[1:], []
        else:
            buf.append(line.strip())
    if header is not None:
        out.append((header, "".join(buf)))
    return out


def header_chain_ids(header):
    """'Chains A, C, L' / 'Chain B[auth H]' / 'Chains E[auth G], F[auth I]' → list of IDs."""
    parts = header.split("|")
    field = (parts[1] if len(parts) > 1 else header).strip()
    field = re.sub(r"^Chains?\s+", "", field, flags=re.I)
    ids = []
    for tok in field.split(","):
        tok = tok.strip()
        if not tok:
            continue
        m = re.search(r"\[auth\s+([A-Za-z0-9]+)\]", tok) or re.match(r"([A-Za-z0-9]+)", tok)
        if m:
            ids.append(m.group(1))
    return ids


def write_fv_fasta(fa_in, out_fa, out_stem, H, L, ag, scheme, rename_map):
    """Cut the RCSB FASTA at the same ANARCI boundary and write it with headers that
    mp_read_fasta_hl.py can parse (HEAVY / LIGHT keywords)."""
    rec_by_chain = {}
    for hdr, seq in parse_fasta(fa_in):
        for cid in header_chain_ids(hdr):
            rec_by_chain.setdefault(cid, (hdr, seq))
    fa_lines, n = [], 0
    for ch, desc in ((H, "Fv heavy chain (VH)"), (L, "Fv light chain (VL)")):
        if ch not in rec_by_chain:
            raise SystemExit(f"[FATAL] no FASTA record for chain '{ch}' "
                             f"(available: {sorted(rec_by_chain)})")
        _, seq = rec_by_chain[ch]
        cut, ctype = anarci_domain(seq, scheme, f"FASTA chain {ch}")
        print(f"  [FASTA] chain {ch}: {len(seq)} → {cut} residues  chain_type={ctype}")
        n += 1
        fa_lines.append(f">{out_stem}_{n}|Chain {rename_map[ch]}|{desc}|MimoPrism preprocess\n")
        fa_lines.append(seq[:cut] + "\n")
    for ch in ag:
        if ch not in rec_by_chain:
            print(f"  [FASTA][warn] no record for antigen chain '{ch}' — skipped", file=sys.stderr)
            continue
        _, seq = rec_by_chain[ch]
        n += 1
        fa_lines.append(f">{out_stem}_{n}|Chain {rename_map[ch]}|antigen (epitope)|MimoPrism preprocess\n")
        fa_lines.append(seq + "\n")
    Path(out_fa).write_text("".join(fa_lines))


# ── 3)+4) Remove water/ligands + renumber (Biopython) ───────────────────────
def is_water(residue):
    return residue.id[0] == "W" or residue.get_resname().strip() in WATER


def is_ligand(residue):
    return residue.id[0] not in (" ", "W") and residue.get_resname().strip() not in WATER


def res_label(chain_id, residue):
    _, resseq, icode = residue.id
    icode = icode.strip()
    return f"{chain_id}{resseq}{icode}" if icode else f"{chain_id}{resseq}"


def res_type(residue):
    if residue.id[0] == " ":
        return "ATOM/protein"
    return "HETATM/water" if is_water(residue) else "HETATM/ligand"


def clean_and_renumber(src_pdb, out_pdb, remove_water, remove_ligand, preserve,
                       renumber_standard_only, write_mapping):
    """src_pdb → remove water/ligands → renumber per chain from 1 → out_pdb (+ mapping CSV).
    Returns the number of remaining protein residues per chain {chain: n}."""
    from Bio.PDB import PDBParser, PDBIO, Select

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("in", str(src_pdb))
    model = structure[0]
    tot = wat = lig = prot = 0
    for chain in model:
        for r in chain:
            tot += 1
            if is_water(r):
                wat += 1
            elif is_ligand(r):
                lig += 1
            else:
                prot += 1
    print(f"  [clean] total {tot} | protein {prot} | water {wat}→{'removed' if remove_water else 'kept'} "
          f"| ligand {lig}→{'removed' if remove_ligand else 'kept'}"
          + (f" | preserved {sorted(preserve)}" if preserve else ""))

    class Clean(Select):
        def accept_residue(self, residue):
            if residue.get_resname().strip().upper() in preserve:
                return True
            if is_water(residue):
                return not remove_water
            if is_ligand(residue):
                return not remove_ligand
            return True

    # Save the cleaned structure and re-parse it so renumbering cannot collide with removed residues.
    tmp = str(Path(out_pdb).with_suffix("")) + "_cleaned_temp.pdb"
    io = PDBIO(); io.set_structure(structure); io.save(tmp, Clean())
    structure = parser.get_structure("clean", tmp)
    model = structure[0]

    rows, counts = [], {}
    for chain in model:
        # Detach every residue first, then re-add with the new id: assigning ids in place
        # collides with the next sibling (e.g. 52A → 53 while 53 still exists) and makes
        # Biopython emit a warning per residue.
        residues = list(chain)
        for r in residues:
            chain.detach_child(r.id)
        n = 1
        for r in residues:
            old = res_label(chain.id, r); rname = r.get_resname(); rtype = res_type(r)
            if (not renumber_standard_only) or r.id[0] == " ":
                r.id = (r.id[0], n, " ")
                rows.append([chain.id, rname, rtype, old, f"{chain.id}{n}", "renumbered"])
                n += 1
            else:
                rows.append([chain.id, rname, rtype, old, old, "preserved"])
            chain.add(r)
        counts[chain.id] = sum(1 for r in chain if r.id[0] == " ")
    io = PDBIO(); io.set_structure(structure); io.save(str(out_pdb))
    os.remove(tmp)
    print(f"  [done] wrote {out_pdb}  ({len(rows)} residues remain)")
    if write_mapping:
        mp = str(Path(out_pdb).with_suffix("")) + "_mapping.csv"
        with open(mp, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["chain", "resname", "residue_type", "original_residue",
                        "renumbered_residue", "status"])
            w.writerows(rows)
        print(f"  [done] mapping: {mp}")
    return counts


# ── main ────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description="raw PDB/FASTA → noncanonical→canonical + Fv cut (ANARCI) + water/ligand removal + renumber")
    ap.add_argument("pdb_name", help="e.g. 6AXK → <dir>/6AXK.pdb, <dir>/rcsb_pdb_6AXK.fasta")
    ap.add_argument("--dir", default=str(DEFAULT_DIR),
                    help="input/output directory (default: <MimoPrism>/input)")
    ap.add_argument("--fasta", default=None, help="input FASTA (default <dir>/rcsb_pdb_<PDB>.fasta)")
    ap.add_argument("--out_dir", default=None, help="output directory (default: --dir)")
    # Fv
    ap.add_argument("--ab", nargs=2, metavar=("CHAIN1", "CHAIN2"), default=None,
                    help="the two antibody chain IDs (heavy/light decided by ANARCI, order free). Required for the Fv cut")
    ap.add_argument("--ag", nargs="+", default=[], help="antigen chain IDs (kept as is)")
    ap.add_argument("--rename", nargs=2, metavar=("H", "L"), default=["H", "L"],
                    help="output heavy/light chain IDs (default H L — always applied)")
    ap.add_argument("--rename_ag", default="P",
                    help="output ID of the first antigen chain (default P — always applied); "
                         "further antigen chains take the following letters (Q, R, ...)")
    ap.add_argument("--no-fv", dest="fv", action="store_false",
                    help="skip the Fv cut and FASTA: canonicalize + clean + renumber only. No suffix, no rename")
    ap.add_argument("--suffix", default="Fv", help="output file suffix (default Fv)")
    ap.add_argument("--scheme", default="imgt", help="ANARCI numbering scheme (imgt/kabat/chothia/martin)")
    # noncanonical residues
    ap.add_argument("--keep_nonstd", action="store_true",
                    help="leave noncanonical residues untouched (default: convert). "
                         "Note: residues left as HETATM are treated as ligands and removed")
    # cleanup / renumber
    ap.add_argument("--keep_water", action="store_true")
    ap.add_argument("--keep_ligand", action="store_true")
    ap.add_argument("--preserve_hetatm", nargs="*", default=[], help="HETATM resnames to keep (e.g. ZN CA)")
    ap.add_argument("--renumber_standard_only", action="store_true",
                    help="renumber ATOM (protein) residues only (default: everything that remains)")
    ap.add_argument("--no_mapping", action="store_true", help="do not write the mapping CSV")
    args = ap.parse_args()

    d = Path(args.dir).expanduser()
    stem = args.pdb_name
    pdb_in = d / f"{stem}.pdb"
    fa_in = Path(args.fasta).expanduser() if args.fasta else d / f"rcsb_pdb_{stem}.fasta"
    out_dir = Path(args.out_dir).expanduser() if args.out_dir else d
    if not pdb_in.exists():
        raise SystemExit(f"[FATAL] input PDB not found: {pdb_in}")
    if args.fv:
        if not args.ab:
            raise SystemExit("[FATAL] the Fv cut needs --ab <chain1> <chain2> (use --no-fv to skip it)")
        if not fa_in.exists():
            raise SystemExit(f"[FATAL] input FASTA not found: {fa_in}  "
                             "(place the RCSB rcsb_pdb_<PDB>.fasta there or pass --fasta)")
        if len(set(args.rename)) != 2:
            raise SystemExit(f"[FATAL] --rename IDs must differ: {args.rename}")
        # antigen output IDs: --rename_ag, then the following letters, skipping the H/L IDs
        ag_new, letter = [], ord(args.rename_ag.upper())
        while len(ag_new) < len(args.ag):
            cid = chr(letter)
            if cid not in args.rename and cid not in ag_new:
                ag_new.append(cid)
            letter += 1
            if letter > ord("Z"):
                raise SystemExit("[FATAL] ran out of chain letters for antigen renaming")
    preserve = {x.strip().upper() for x in args.preserve_hetatm}

    # Do not append the suffix twice when the name already ends with it (re-runs / manual Fv input).
    if args.fv:
        out_stem = stem if stem.endswith(f"_{args.suffix}") else f"{stem}_{args.suffix}"
    else:
        out_stem = stem

    print(f"[preprocess] input PDB   : {pdb_in}")
    if args.fv:
        print(f"[preprocess] input FASTA : {fa_in}")
    lines = pdb_in.read_text().splitlines(keepends=True)

    # ── 1) noncanonical residues ────────────────────────────────────────────
    unknown = unknown_hetatm(lines)
    if args.keep_nonstd:
        print("[1/4] noncanonical residues: --keep_nonstd — not converted")
    else:
        lines, report, dropped = convert_nonstd(lines)
        lines = clean_records(lines, dropped, {k.split("(")[1].split("→")[0] for k in report})
        if report:
            print("[1/4] noncanonical → canonical residues (HETATM → ATOM):")
            for key, r in report.items():
                chain_res, conv = key.split("(")
                parent = conv.split("→")[1].rstrip(")")
                parts = [f"hydrogens removed: {r['hydrogens']}"]
                if r["extra"]:
                    parts.append(f"extra heavy atoms removed: {','.join(r['extra'])}")
                if r["missing"]:
                    parts.append(f"{parent} atoms absent in the original (not rebuilt): {','.join(r['missing'])}")
                print(f"  {chain_res:<6} {conv.rstrip(')'):<10} {' | '.join(parts)}")
        else:
            print("[1/4] noncanonical → canonical residues: nothing to convert (canonical residues only)")
    if unknown:
        print(f"  [warn] HETATM residues not in the table — treated as ligands: {', '.join(unknown)}", file=sys.stderr)
        print("         (if you don't want, keep a ligand with --preserve_hetatm)", file=sys.stderr)

    out_dir.mkdir(parents=True, exist_ok=True)

    # ── 2) Fv cut (+ FASTA) ─────────────────────────────────────────────────
    H = L = None
    if args.fv:
        print(f"[2/4] Fv cut (ANARCI {args.scheme})")
        lines, H, L, rename_map = cut_fv(lines, args.ab, args.ag, args.scheme, args.rename, ag_new)
        hid, lid = args.rename
        print(f"  [rename] heavy: {H} → {hid},  light: {L} → {lid}"
              + (",  antigen: " + ", ".join(f"{a} → {b}" for a, b in zip(args.ag, ag_new))
                 if args.ag else "  (no antigen chain)"))
        fv_pdb = out_dir / f"{out_stem}.pdb"
        fv_pdb.write_text("".join(lines))
        print(f"  [wrote] {fv_pdb}  (original numbering kept)")
        out_fa = out_dir / f"rcsb_pdb_{out_stem}.fasta"
        write_fv_fasta(fa_in, out_fa, out_stem, H, L, args.ag, args.scheme, rename_map)
        print(f"  [wrote] {out_fa}")
        src_for_clean = fv_pdb
    else:
        print("[2/4] Fv cut: --no-fv — skipped")
        if args.ab or args.ag:                   # chain selection only
            keep = set(args.ab or []) | set(args.ag)
            lines = [l for l in lines
                     if not l.startswith(("ATOM", "HETATM", "ANISOU")) or l[21] in keep]
        src_for_clean = out_dir / f"{out_stem}_std_temp.pdb"
        src_for_clean.write_text("".join(lines))

    # ── 3)+4) remove water/ligands + renumber ───────────────────────────────
    print("[3/4] remove water/ligands  [4/4] renumber per chain")
    processed = out_dir / f"processed_{out_stem}.pdb"
    counts = clean_and_renumber(src_for_clean, processed,
                                remove_water=not args.keep_water,
                                remove_ligand=not args.keep_ligand,
                                preserve=preserve,
                                renumber_standard_only=args.renumber_standard_only,
                                write_mapping=not args.no_mapping)
    if not args.fv:
        os.remove(src_for_clean)

    # ── report ──────────────────────────────────────────────────────────────
    print(f"\n[preprocess] output name: {out_stem}   (use as PDB_NAME)")
    if args.fv:
        hid, lid = args.rename
        nh, nl = counts.get(hid, 0), counts.get(lid, 0)
        print(f"[preprocess] receptor: VH={nh}({hid}) + VL={nl}({lid}) residues (processed PDB)")
        print("\nnext step:")
        print(f"  python preprocess/select_hotspot.py {out_stem} --ab {hid} {lid}"
              + (f" --ag {' '.join(ag_new)}" if ag_new else "") + " --ks 6")
        print(f"  CONTIGS_RECEPTOR=\"{hid}1-{nh}/0 {lid}1-{nl}/0\"")


if __name__ == "__main__":
    main()