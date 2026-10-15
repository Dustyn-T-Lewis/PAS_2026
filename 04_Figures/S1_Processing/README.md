# S1_Processing

Supplementary Figure 1: what processing kept and chose.

![Supplementary Figure 1](b_reports/S1_composite.png)

| | |
|---|---|
| Reads | `01_Processing/01_Filtering/c_data/01_filtering.xlsx`, `01_Processing/02_Quantification/c_data/02_quantification.xlsx`, `01_Processing/03_Design/c_data/03_design.xlsx` |
| Writes | `b_reports/panels/<panel>.pdf/.png`, `b_reports/S1_composite.pdf/.png` (180 × 117 mm), `c_data/S1.xlsx` |
| Run | `Rscript 04_Figures/S1_Processing/a_script/S1_composite.R` (each panel script also runs alone) |

Each panel is its own script in `a_script/panels/`; the composite runs them, adds the letters, stitches them and checks the mouse counts and that each chosen normalisation was eligible. Panels read saved outputs and compute nothing beyond counts, filtering and scaling.

## Panels

| Panel | Shows | Source |
|---|---|---|
| a | Precursors and protein groups left after each filter. | `filter_log` |
| b | Precursors and proteins identified per run, by group. | `run_qc` |
| c | Label-blind normalisation benchmark: within-group PCV, PMAD and PEV per method and universe; only methods that keep limpa's standard errors were eligible. | `normalization` |
| d | The surrogate variable against run depth: r 0.85 overall, 0.86 in females, 0.92 in males; no treatment difference. | `surrogate`, `surrogate_r` |
| e | Share of pellet signal from MitoCarta proteins: differs by sex (p 2e-8), not by treatment (p 0.83). | `mito_share`, `mito_share_tests` |

## Workbook

`c_data/S1.xlsx` holds one sheet per panel with the numbers it draws, plus `read_me`, `input_manifest` and `package_versions`.
