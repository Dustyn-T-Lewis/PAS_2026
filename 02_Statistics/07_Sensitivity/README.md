# 07_Sensitivity

How much the protein results depend on choices made once. Each row changes one choice from the baseline and refits the three comparison contrasts in the full proteome exactly as 01_Fit does.

| | |
|---|---|
| Reads | `01_Processing/01_Filtering/c_data/precursors_filtered.rds`, `01_Processing/02_Quantification/c_data/proteins.rds`, `01_Processing/03_Design/c_data/design.rds`, `02_Statistics/01_Fit/c_data/fit.rds` |
| Writes | `c_data/07_sensitivity.xlsx`, `c_data/sensitivity.rds` (read by `04_Figures/F02_Proteins` for S4), `c_data/dpcquant_slope_*.rds` (reusable quantification at other slopes), `b_reports/07_sensitivity.pdf` |
| Run | `Rscript 02_Statistics/07_Sensitivity/a_script/07_sensitivity.R` (about six minutes the first time) |

| Variant | Change from baseline |
|---|---|
| slope 0.5, slope 1.0 | limpa detection-curve slope instead of 0.7 (dpcQuant rerun; same 1,962 proteins; surrogate variable kept from the baseline) |
| run depth instead of sv1 | measured log2 run depth as the covariate |
| no adjustment | no covariate |
| no sample weights | `dpcDE(sample.weights = FALSE)` |
| without PAS_M1 | the re-injected male PAS run removed |
| pellet share covariate | pellet MitoCarta share added to the design |

- **Reported:** per variant and contrast, proteins at p ≤ 0.05, FDR ≤ 0.10 and FDR ≤ 0.05, and how many of the baseline's 50 smallest p-values stay among the variant's 50 smallest; for every protein at FDR ≤ 0.10 in the baseline, its p and FDR under each variant. The baseline row is checked against 01_Fit.
- **Result:** the female top three (Dusp3, Cops3, Cops8) stay at FDR ≤ 0.10 under every variant. All three stay at FDR ≤ 0.05 under run depth, no sample weights, without PAS_M1 and pellet share; slope 0.5 keeps only Cops3 there, slope 1.0 keeps Cops8 and Dusp3, and the unadjusted fit keeps only Dusp3. The four male proteins at FDR ≤ 0.10 in the baseline (Cd36, Cd34, Nufip2, Slc2a4) fall above it under slope 0.5, no adjustment, no sample weights and without PAS_M1. No variant gives an interaction protein at FDR ≤ 0.10. Counts at p ≤ 0.05 move little except in males without adjustment or weights (193 to 135 and 132), where run depth goes unhandled. Without sample weights, females reach 13 proteins at FDR ≤ 0.05.
