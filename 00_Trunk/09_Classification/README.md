# 09_Classification

Asks how well each pathway score, module eigengene, complex score and cluster score alone separates PAS from vehicle mice within each sex.

- **Method:** `pROC::roc` AUC (1 = higher in PAS) with DeLong 95% CI, left empty at perfect separation; exact Wilcoxon p; BH within sex and level; r with group beside each row. The number of perfect separators is compared with all 252 relabellings of the same ten mice.
- **Inputs:** `07_Fit/c_data/fit.rds`
- **Outputs:** `c_data/09_classification.xlsx`, `c_data/classification.rds` (read by the female and male `04_Classification`), `b_reports/09_classification.pdf`
- **Result:** female: 47 of 889 pathway scores at p ≤ 0.05 (44 expected), 17 perfect separators (null median 3, p = 0.095). Male: 20 of 889 (44 expected), 3 perfect (p = 0.41). No feature reaches FDR ≤ 0.05 in either sex.
- **Note:** proteins are left out because their ROC repeats the protein contrast. The sex-difference arc has no classification; its test is the interaction contrast.
