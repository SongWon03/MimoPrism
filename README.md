# MimoPrism (v0.5)

## What is MimoPrism?  

A mimotope is a peptide that mimics the conformation of an epitope-the specific region of a protein that interacts with antibodies or T-cell receptors. We are focusing on designing mimotopes that can bind to the complementarity-determining regions (CDRs) of [monoclonal antibodies (mAbs)](https://en.wikipedia.org/wiki/Monoclonal_antibody) like natural epitope. These novel antigens may elicit high-affinity antibodies upon immunization, potentially reducing costs and enhancing efficacy compared to the passive administration of expensive mAbs.  

**MimoPrism** is a computational pipeline for mimotope design and *in silico* evaluation. The input of MimoPrism is the experimentally determined complex structure of peptide antigen and mAb. And the output is a list of designed mimotopes with their predicted complex structure and scores.  


## MimoPrism Pipeline  

You can run our pipeline at once with **run_mimoprism.sh**.  
The tunable parameters for pipeline can be selected by **mimoprism_config.json**.  

**`RFpeptides` &rarr; `CyclicMPNN` &rarr; `ESMFold2 (no MSA)` &rarr; `AlphaFold3 (MSA)`**  

1. Generate N cyclic peptide [backbone](https://proteopedia.org/w/Backbone_representations) designs for each parameter set (**hotspot residues**, **length**) with `RFpeptide`. 
2. Generate M sequence designs (random temperature in range 0.1-0.3, [recommended](https://docs.nvidia.com/nim/bionemo/proteinmpnn/latest/endpoints.html)) for each backbone by `CycliMPNN`.
3. Based on `CyclicMPNN`'s global score (lower is best), select best sequence for each backbone. Now we have N X M designed sequences for each parameter set. 
4. Predict the complex structure of target mAb and designed mimotope candidates with `ESMFold2 (no MSA)`. 
5. Select top K structures based on `ESMFold2 (no MSA)` **ipTM** score for each parameter set, and cross-validate the complex structure with `AlphaFold3`.
6. Combine and filter the results based on our scoring functions.


## Hotspot Residues Selection

We use two strategies for hotspot selection. 

1. **Beta carbon distance between natural antigen-monoclonal antibody complex**: From the [RFdiffusion Github document](https://github.com/RosettaCommons/RFdiffusion?tab=readme-ov-file#macrocyclic-peptide-design-with-rfpeptides), the definition of **hotspot residues** is part of residues in binding interface that binder must contact with target protein. In training step, they selected all residues which are within 10Å Cbeta distance of binder and randomly used 0-20% of them as interface hotspots. Based on this, we first compute the all residues in mAb which are within 10Å Cbeta distance of original peptide antigen and select nearest 3-7 residues for hotspot.
2. **Paratope prediction by ParaSurf**: [ParaSurf](https://github.com/aggelos-michael-papadopoulos/ParaSurf) is deep learning based program to predict the paratope of given antibody.


## Scoring Functions  

**Metrics**  

- **ESMFold2 (no MSA)/ipTM**: ipTM (interface predicted TM-score) measures the accuracy of the predicted relative positions of the the subunits in a complex structure. Values higher than 0.8 represent confident high-quality predictions, while values below 0.6 suggest a failed prediction. ipTM values between 0.6 and 0.8 are a gray zone where predictions could be correct or incorrect. MimoCraft utilize ipTM as a primary filtering criterion for evaluating confidence of the predicted complex structure.  
- **ESMFold2 (no MSA)/pLDDT**: pLDDT (position local distance difference test) is confidence metric of the position of amino acid residues. > 0.9 means good.  
- **AlphaFold3/ranking_score**: A scalar in the range [-100, 1.5] that corporates `ptm`, `iptm`, `fraction_disordered`, and `has_clash` into a single number with the following equation: 0.8 × ipTM + 0.2 × pTM + 0.5 × disorder − 100 × has_clash. Higher is better.
- **AlphaFold3/ipTM**: ipTM (interface predicted TM-score) > 0.8 means good.  
- **AlphaFold3/pLDDT**: pLDDT (position local distance difference test) is confidence metric of the position of amino acid residues. > 0.9 means good.  
- **AlphaFold3/iPAE**: iPAE (interface predicted aligned error) < 0.4 (generally used)
- **Rosetta**: Calculates **binding affinity**, **spatial aggregation propensity**, **contact molecular surface area (CMS)**.
