# 03_Classify

Asks how well each protein on its own separates the two groups of each comparison. Descriptive:
with five mice per group, perfect separation happens by chance with probability 2/252.

- **Method:** `pROC::roc()` AUC with DeLong 95% interval (`ci.auc(method = "delong")`, undefined at perfect separation) and exact `wilcox.test()`, BH within comparison. Four comparisons: `PAS_vs_VEH_female`, `PAS_vs_VEH_male`, `sex_in_VEH`, `sex_in_PAS` (first group minus second). No interaction comparison.
- **Inputs:** `02_Differential_Expression/02_Contrasts/c_data/contrasts.rds` (quantile-normalised full universe, 1,842 proteins)
- **Outputs:** `c_data/03_classify.xlsx` (sheets: summary, auc, plot_index), `b_reports/03_classify.pdf` (strip charts of the 12 strongest separators, one page per comparison)
- **Result:** 14.6 perfect separations expected by chance per comparison. Female PAS: 41 perfect, 117 at p ≤ 0.05, none at FDR ≤ 0.10. Male PAS: 7 perfect, 52 at p ≤ 0.05, none at FDR ≤ 0.10. Sex in VEH: 305 perfect, 414 at FDR ≤ 0.10. Sex in PAS: 492 perfect, 743 at FDR ≤ 0.10.
- **Note:** The smallest exact Wilcoxon p at 5 against 5 is 0.008, so no treatment comparison can reach FDR ≤ 0.10 across 1,842 proteins unless many proteins separate perfectly.
