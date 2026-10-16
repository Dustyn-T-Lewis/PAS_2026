# S8_Sensitivity

Supplementary Figure 8: how the protein results move when one processing or model choice changes.

![Supplementary Figure 8](b_reports/S8_composite.png)

| | |
|---|---|
| Reads | `02_Statistics/07_Sensitivity/c_data/sensitivity.rds`, `07_sensitivity.xlsx` |
| Writes | `b_reports/panels/<panel>.pdf/.png`, `b_reports/S8_composite.pdf/.png` (180 × 162 mm), `c_data/S8.xlsx` |
| Run | `Rscript 04_Figures/S8_Sensitivity/a_script/S8_composite.R` (each panel script also runs alone) |

Each variant changes one choice and refits the three comparison contrasts in the full proteome (see `02_Statistics/07_Sensitivity/README.md`).

## Panels

| Panel | Shows | Source |
|---|---|---|
| a | Counts per variant and contrast at p ≤ 0.05, FDR ≤ 0.10 and FDR ≤ 0.05, and how many of the baseline's 50 smallest p-values stay in the variant's 50. Cell colour: log2 ratio to the baseline, one added to both counts. Every variant finds 132 to 212 proteins at p ≤ 0.05 where chance alone gives 98. The female FDR ≤ 0.10 count runs from 4 (no adjustment) to 25 (run depth instead of sv1). The male FDR ≤ 0.10 count is 4 at baseline and 0 under four variants. No variant finds an interaction protein at FDR ≤ 0.10. | `counts` |
| b | The 22 proteins at FDR ≤ 0.10 in the baseline, binned FDR under each variant, with the log2 FC as text. Cops3, Dusp3 and Cops8 stay at FDR ≤ 0.10 under every variant and never change sign. Of the four male proteins, Nufip2 and Cd34 stay at FDR ≤ 0.10 under 3 of the 7 variants, Cd36 under 1 and Slc2a4 under none. | `named` |

## Workbook

`c_data/S8.xlsx` holds one sheet per panel with the numbers it draws, plus `read_me`, `input_manifest` and `package_versions`.
