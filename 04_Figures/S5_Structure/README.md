# S5_Structure

Supplementary Figure 5: whether perfect separation exceeds chance, and what the proteome components track.

![Supplementary Figure 5](b_reports/S5_composite.png)

| | |
|---|---|
| Reads | `02_Statistics/03_Classification/c_data/classification.rds`, `02_Statistics/04_PCA/c_data/pca.rds` |
| Writes | `b_reports/panels/<panel>.pdf/.png`, `b_reports/S5_composite.pdf/.png` (180 × 60 mm), `c_data/S5.xlsx` |
| Run | `Rscript 04_Figures/S5_Structure/a_script/S5_composite.R` (each panel script also runs alone) |

Each panel is its own script in `a_script/panels/`; the composite runs them, adds the letters, stitches them and checks the observed counts against `03_classification.xlsx`. Panels read saved outputs and compute nothing beyond counts, filtering and scaling.

## Panels

| Panel | Shows | Source |
|---|---|---|
| a | Features that put all five PAS mice above or below all five vehicle mice, per level and sex, against all 252 relabellings: only female modules exceed chance (2 against a null median of 0, p 0.016). | `classification.rds`: `perfect` |
| b | R² of sex, treatment, their interaction, MitoCarta share and run depth on proteome PC1 to PC4: PC1 tracks sex and MitoCarta share; treatment reaches 0.21 only on PC4. | `pca.rds`: `component_r2` |

## Workbook

`c_data/S5.xlsx` holds one sheet per panel with the numbers it draws, plus `read_me`, `input_manifest` and `package_versions`.
