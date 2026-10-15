# F02_Proteins

Figure 2: which proteins PAS changes in either sex, and how the two sexes' responses line up.

![Figure 2](b_reports/F02_composite.png)

| | |
|---|---|
| Reads | `02_Statistics/01_Fit/c_data/fit.rds`, `02_Statistics/02_Enrichment/c_data/enrichment.rds` |
| Writes | `b_reports/panels/<panel>.pdf/.png`, `b_reports/F02_composite.pdf/.png` (180 × 145 mm), `c_data/F02.xlsx` |
| Run | `Rscript 04_Figures/F02_Proteins/a_script/F02_composite.R` (each panel script also runs alone) |

Each panel is its own script in `a_script/panels/`; the composite runs them, adds the letters, stitches them and checks the heatmap rows and pattern counts against `01_fit.xlsx` and the ORA FDRs against `02_enrichment.xlsx`. Panels read saved outputs and compute nothing beyond counts, filtering and scaling.

## Panels

| Panel | Shows | Source |
|---|---|---|
| a | The 215 full-proteome proteins at moderated F p ≤ 0.05 over both PAS contrasts (none at FDR ≤ 0.05), across all 20 mice, each centred on its own sex's vehicle mean; rows split by response pattern (56 up in both, 83 down in both, 56 female up / male down, 20 female down / male up) and clustered within each split. Right: each protein's pooled r with each phenotype within each sex. The rows were chosen for differing, so the clustering describes them and tests nothing. | `fit.rds`: `protein_ftest`, `E_norm`, `protein_associations` |
| b | Over-representation (fgsea::fora) of each pattern against the detected proteome, full collections: the three smallest-p sets per pattern, outlined at FDR ≤ 0.05. Only the female up / male down pattern has sets at FDR ≤ 0.05 (Reactome striated muscle contraction and muscle contraction). | `enrichment.rds`: `pattern_ora` |

## Workbook

`c_data/F02.xlsx` holds one sheet per panel with the numbers it draws, plus `read_me`, `input_manifest` and `package_versions`.
