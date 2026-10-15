# 03_Classification

Asks how well each pathway score, module eigengene, complex score and cluster score alone separates PAS from vehicle mice within each sex.

| | |
|---|---|
| Reads | `02_Statistics/01_Fit/c_data/fit.rds` |
| Writes | `c_data/03_classification.xlsx`, `c_data/classification.rds` (read by the female and male `04_Classification` and `04_Figures/`), `b_reports/03_classification.pdf` |
| Run | `Rscript 02_Statistics/03_Classification/a_script/03_classification.R` |

- **Method:** `pROC::roc` AUC (1 = higher in PAS) with DeLong 95% CI, left empty at perfect separation; exact Wilcoxon p; BH within sex and level; r with group beside each row. The number of perfect separators is compared with all 252 relabellings of the same ten mice. Each module's ROC curve (`pROC::coords`, sheet `module_curves`) is saved for the figures.
- **Result:** female: 47 of 1,005 pathway scores at p ≤ 0.05 (50 expected), 16 perfect separators (null median 4, p 0.13); 2 of 14 modules separate perfectly (null median 0, p 0.016). Male: 28 of 1,005 (50 expected), 5 perfect (null median 2, p 0.33). No feature reaches FDR ≤ 0.05 in either sex.
- **Basis:** Robin et al. 2011 (pROC) for AUC and the DeLong interval; the AUC equals the Mann-Whitney U over n₁n₀, and the perfect-separator null enumerates all choose(10, 5) = 252 labellings, both recomputed by hand to match. Scores enter unadjusted: a single-feature AUC has no model to carry the surrogate variable.
- **Note:** Proteins are left out because their ROC repeats the protein contrast. The sex-difference comparison has no classification; its test is the interaction contrast.
