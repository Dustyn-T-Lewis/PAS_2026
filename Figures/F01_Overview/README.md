# F01_Overview

Figure 1: the design and the three questions side by side (PAS vs vehicle in females, in males, and the treatment × sex interaction) before any single result.

| | |
|---|---|
| Reads | `assets/design.png`, `00_Trunk/07_Fit/c_data/fit.rds`, `00_Trunk/08_Enrichment/c_data/enrichment.rds`, `00_Trunk/10_Global/c_data/global.rds` |
| Writes | `b_reports/f01_overview.pdf` and `.png` (300 × 190 mm), `c_data/f01_overview.xlsx` |
| Run | `Rscript Figures/F01_Overview/a_script/f01_overview.R` |

The script reads saved outputs only and computes nothing beyond counts and arrow scaling. Before writing, it checks the panel b FDR counts against `protein_contrasts`, the mouse counts, that every pathway row has all three columns, and that every selected set is at FDR ≤ 0.05 in its source.

## Panels

| Panel | Shows | Source |
|---|---|---|
| a | Design: 2 × 2 sex by treatment, n = 5 per cell; gastrocnemius mitochondrial pellet, DIA; the phenotypes; the three contrasts in the colours used throughout. | `assets/design.png` |
| b | Proteins per contrast, down left and up right, full proteome and mitochondrial universe. Light bar: p ≤ 0.05 or Π ≤ 0.05; dark bar, drawn over it from zero: Π ≤ 0.05. Dashed lines: p ≤ 0.05 expected by chance on each side. FDR ≤ 0.05 count beside each row. | `fit.rds`: `protein_contrasts` |
| c | Proteome PCA of the 20 mice, a 95% t-ellipse per group, PERMANOVA, and the MitoCarta share fitted as a vector (R² 0.89). | `global.rds` |
| d | Phenotype PCA of the six traits, scaled, with loadings and PERMANOVA. | `global.rds` |
| e | Phenotype effects: PAS vs vehicle within each sex, % change with 95% CI on a log axis (Welch t on log values); the interaction p under each trait. | `fit.rds`: `phenotype_effects`, `phenotype_interaction` |
| f | Every set at fgsea FDR ≤ 0.05 in any contrast that the EnrichmentMap flag keeps, listed once under the first contrast that selects it and drawn in all three. Square: collection. † : the same set's singscore contrast, tested across mice, has p > 0.05. | `enrichment.rds`: `results` |
| g | Each sex's p ≤ 0.05 up and down proteins placed in the other sex's ranking, as limma's barcode (Wu & Smyth 2012); each lane gives size, the share ranked the same way and the fgsea p. | `enrichment.rds`: `signatures`, `ranks` |

Π is p^|log2 FC| (Xiao et al. 2014). It ranks a protein on p and fold change together and controls no error rate; Π ≤ 0.05 can include proteins with p > 0.05.

Panel a was generated with Nano Banana from `assets/design_spec/` (SPEC, EVIDENCE, DONOTDRAW; prompt and model in PROMPT.txt and PROVENANCE.json) and checked label by label. Known styling defects were accepted: the column connectors cross the processing row, the male-to-interaction connector is blue, and the purple connector is doubled. To change it, edit the spec and rerun `build_from_spec.py`, or replace `design.png` with a BioRender export; nothing else changes.

## Workbook

| Sheet | Holds |
|---|---|
| `protein_counts` | Panel b counts per universe, contrast and direction. |
| `pca_proteome`, `pca_phenotype` | Panels c, d scores per mouse. |
| `permanova` | Panels c, d tests. |
| `phenotype_effects`, `phenotype_interaction` | Panel e. |
| `pathways` | Panel f: each selected set in all three contrasts, with NES, p, FDR, singscore p and r with group. |
| `signatures` | Panel g lanes. |
