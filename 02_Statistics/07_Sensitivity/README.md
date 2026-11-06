# 07_Sensitivity

Tests how much the protein results depend on choices made once. Each variant changes one choice and refits the three comparison contrasts in the full proteome as 01_Fit does.

| | |
|---|---|
| Reads | `01_Processing/01_Filtering/c_data/precursors_filtered.rds`, `01_Processing/02_Quantification/c_data/proteins.rds`, `01_Processing/03_Design/c_data/design.rds`, `02_Statistics/01_Fit/c_data/fit.rds` |
| Writes | `c_data/07_sensitivity.xlsx`, `c_data/sensitivity.rds`, `c_data/dpcquant_slope_0.5.rds`, `c_data/dpcquant_slope_1.rds`, `b_reports/07_sensitivity.pdf` |
| Run | `Rscript 02_Statistics/07_Sensitivity/a_script/07_sensitivity.R` |

## Method

| Variant | Change from baseline |
|---|---|
| slope 0.5, slope 1.0 | limpa detection-curve slope instead of 0.7 (same 1,962 proteins; baseline surrogate variable kept) |
| run depth instead of sv1 | measured log2 run depth as the covariate |
| no adjustment | no covariate |
| no sample weights | `dpcDE(sample.weights = FALSE)` |
| without PAS_M1 | the re-injected male PAS run removed |
| pellet share covariate | pellet MitoCarta share added to the design |

- Per variant and contrast: proteins at p ≤ 0.05, FDR ≤ 0.10 and FDR ≤ 0.05, and how many of the baseline's 50 smallest p stay among the variant's 50 smallest. Every baseline protein at FDR ≤ 0.10 is followed across variants.

## Results

- Cops3, Dusp3 and Cops8 stay at FDR ≤ 0.10 under every variant, and at FDR ≤ 0.05 under five, six and five of the seven.
- None of the four male proteins at baseline FDR ≤ 0.10 (Cd34, Cd36, Nufip2, Slc2a4) holds under every variant; Slc2a4 holds under none.
- No variant gives an interaction protein at FDR ≤ 0.10. No baseline protein changes sign.
- Counts at p ≤ 0.05 move little except in males without adjustment (135) or sample weights (132), against 193 at baseline.
