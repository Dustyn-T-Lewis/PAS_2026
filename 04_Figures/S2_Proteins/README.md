# S2_Proteins

Supplementary Figure 2: each comparison's own protein list.

![Supplementary Figure 2](b_reports/S2_composite.png)

| | |
|---|---|
| Reads | `02_Statistics/01_Fit/c_data/fit.rds` |
| Writes | `b_reports/panels/<panel>.pdf/.png`, `b_reports/S2_composite.pdf/.png` (180 × 192 mm), `c_data/S2.xlsx` |
| Run | `Rscript 04_Figures/S2_Proteins/a_script/S2_composite.R` (each panel script also runs alone) |

Each panel is its own script in `a_script/panels/`; the composite runs them, adds the letters, stitches them and checks each heatmap's row count against the Π counts in `01_fit.xlsx`. Panels read saved outputs and compute nothing beyond counts, filtering and scaling.

## Panels

| Panel | Shows | Source |
|---|---|---|
| a | Female PAS vs vehicle: the 107 proteins at Π ≤ 0.05 across all mice, centred on own-sex vehicle, split by direction. | `fit.rds`: `protein_contrasts`, `E_norm` |
| b | Male PAS vs vehicle: 106 proteins at Π ≤ 0.05. | as a |
| c | Treatment × sex: 142 proteins at Π ≤ 0.05, split by which sex's PAS effect is larger. | as a |
| d | Every full-proteome protein per contrast: log2 FC against −log10 p, Π ≤ 0.05 coloured by direction, FDR ≤ 0.05 ringed (three, all female). | `fit.rds`: `protein_contrasts` |

Π is p^|log2 FC| (Xiao et al. 2014); it ranks on p and fold change together and controls no error rate.

## Workbook

`c_data/S2.xlsx` holds one sheet per panel with the numbers it draws, plus `read_me`, `input_manifest` and `package_versions`.
