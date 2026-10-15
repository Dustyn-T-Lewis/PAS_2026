# S4_Protein_associations

Supplementary Figure 4: how many proteins track each phenotype, and the strongest one.

![Supplementary Figure 4](b_reports/S4_composite.png)

| | |
|---|---|
| Reads | `02_Statistics/01_Fit/c_data/01_fit.xlsx`, `02_Statistics/01_Fit/c_data/fit.rds` |
| Writes | `b_reports/panels/<panel>.pdf/.png`, `b_reports/S4_composite.pdf/.png` (180 × 135 mm), `c_data/S4.xlsx` |
| Run | `Rscript 04_Figures/S4_Protein_associations/a_script/S4_composite.R` (each panel script also runs alone) |

Each panel is its own script in `a_script/panels/`; the composite runs them, adds the letters, stitches them and checks the counts and that each shown protein appears under both models. Panels read saved outputs and compute nothing beyond counts, filtering and scaling.

## Panels

| Panel | Shows | Source |
|---|---|---|
| a | Proteins at p ≤ 0.05 per phenotype and sex under the pooled, pooled + MitoCarta-share and within-treatment models, against the 98 expected by chance. | `01_fit.xlsx`: `association_summary` |
| b | For each phenotype and sex, the protein with the smallest within-treatment p, and its r-equivalent under the pooled and within-treatment models; filled at FDR ≤ 0.05. | `fit.rds`: `protein_associations` |

## Workbook

`c_data/S4.xlsx` holds one sheet per panel with the numbers it draws, plus `read_me`, `input_manifest` and `package_versions`.
