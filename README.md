# PAS_2026

Gastrocnemius proteomics of pipecolic acid plus succinate (PAS) against vehicle in mice with a
ligated femoral artery. 20 mice, 5 per sex and treatment, one mitochondria-enriched pellet each,
DIA mass spectrometry, 6 muscle phenotypes. 1,842 proteins, 568 in MitoCarta.

Contrasts: PAS vs vehicle in females, in males, their interaction, and sex within each treatment.

## Results

- Sex drives the proteome: PERMANOVA R² 0.38, p = 0.001. Treatment: R² 0.04, p = 0.27.
- PAS moves 152 female and 133 male proteins at p ≤ 0.05 (92 expected by chance). FDR ≤ 0.05:
  4 in females, none in males.
- Females: muscle contraction, myogenesis and myofibril pathways up; cell cycle down.
  Males: OXPHOS, EMT and actin pathways down.
- The two responses are independent (fry p > 0.5). Contractile proteins rise with PAS in females
  and fall in males.
- PAS raises capillary density, CSA and satellite cells in females. The protein signature does
  not track these traits within groups.

## Layout

```
00_Input/                    data
01_Preprocess/               filtering, limpa quantification
02_Differential_Expression/  design, contrasts, classification, phenotype association
03_Pathway_Enrichment/       gene sets, scores, contrasts, classification, association, plots,
                             compartments, figure term selection
04_Network/                  WGCNA modules, contrasts, classification, association, STRING, INDRA
05_Summary/                  signal, panels, biclusters, sex patterns, global structure,
                             magnitude, signature carry-over, signature against phenotype
06_Figures/                  F01-F05, supplements S1-S4
R/                           helpers.R, figures.R, panels.R
Supplementary/               submission files, once drafted
run_all.R                    every step in order
```

Each sub-stage: `a_script/` (one script), `b_reports/` (PDFs), `c_data/` (one workbook, plus
`.rds` read downstream). Stages 01 to 05 are exploratory. Workbooks give p, FDR and Π side by side
with a contrast column on every row, open with `read_me` and end with input checksums and package
versions.

Figures: `a_script/main/F0x.R` stitches the panel scripts in `main/panels/`; supplements sit in
`supp/` the same way; `b_reports/` mirrors this; `c_data/` holds one workbook per figure.

## Rules

- Significant: p ≤ 0.05 or FDR ≤ 0.05 (BH within contrast and collection). Proteins also at
  FDR ≤ 0.10. Π = p^|log2 FC| ranks; it is not an error rate.
- A sex difference needs the interaction contrast. A hit in one sex alone is "detected in" that sex.
- Phenotypes come from soleus, plantaris and IHC, the proteome from gastrocnemius: associations only.

## Run

```sh
Rscript -e 'renv::restore()'
Rscript run_all.R
```

R 4.6.0; `renv.lock` pins every package, enrichVolcano from
github.com/Dustyn-T-Lewis/enrichVolcano. STRING, INDRA and the classification panels cache their
results; rebuilding the panel cache takes about 25 minutes on 12 cores.

## Reference

Burke BI, Valentino TR, Ismaeel A, et al. Exercise-associated microbial metabolites prevent
skeletal muscle atrophy in adult female mice. Nat Commun 2026. PMID 42431880.
