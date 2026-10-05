# 03_Classify

How well each set score alone separates the two groups of each comparison. Descriptive: it ranks
sets, it does not test pathways.

- **Method:** `auc_screen()` in `R/helpers.R`: `pROC` AUC with DeLong 95% interval, exact Wilcoxon
  p, BH within comparison. Four comparisons (female PAS, male PAS, sex in VEH, sex in PAS; no
  interaction).
- **Inputs:** `01_Scores/c_data/set_scores.rds`,
  `02_Differential_Expression/02_Contrasts/c_data/contrasts.rds`
- **Outputs:** `c_data/03_classify.xlsx` (sheets: summary, auc, plot_index),
  `b_reports/03_classify.pdf` (12 strongest separators per comparison)
- **Result:** 949 sets per comparison; 7.5 perfect separations (AUC 0 or 1) expected by chance.
  Observed: female PAS 18, male PAS 4, sex in VEH 103, sex in PAS 186. Only sex in PAS has sets at
  FDR ≤ 0.05 (the 186 perfect separators); treatment comparisons have none below FDR 0.10.
- **Note:** With 5 against 5 the smallest exact Wilcoxon p is 2/252 = 0.0079, so FDR ≤ 0.05 is
  reachable only when many sets separate perfectly.
