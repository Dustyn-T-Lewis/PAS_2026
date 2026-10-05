# 03_Classify

Asks how well each module eigengene alone separates the two groups of each comparison. Descriptive.

- **Method:** shared `auc_screen()` in `R/helpers.R`: `pROC` AUC with DeLong 95% interval, exact Wilcoxon p, BH within comparison; four comparisons (PAS vs VEH in each sex, sex in VEH, sex in PAS).
- **Inputs:** `../01_Build_Modules/c_data/modules.rds`, `02_Differential_Expression/02_Contrasts/c_data/contrasts.rds`
- **Outputs:** `c_data/03_classify.xlsx` (sheets: summary, auc, plot_index), `b_reports/03_classify.pdf`
- **Result:** no module separates PAS from VEH in females (0 of 9 at p ≤ 0.05). In males only black does (AUC 0.08, p = 0.032, FDR 0.29). By sex, turquoise and blue separate perfectly in VEH (p = 0.008, FDR 0.036), and turquoise, blue, brown and pink in PAS (FDR 0.018).
- **Note:** at 5 against 5 the smallest exact p is 0.0079 and chance gives a perfect split with probability 2/252; DeLong intervals collapse at perfect separation.
