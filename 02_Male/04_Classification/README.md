# 04_Classification

Does one pathway score, module eigengene, complex score or cluster score alone separate PAS from vehicle mice in males, five against five?

| | |
|---|---|
| Reads | `00_Trunk/09_Classification/c_data/classification.rds` |
| Writes | `c_data/04_classification.xlsx`, `c_data/classification.rds`, `b_reports/04_classification.pdf` |
| Run | `Rscript 02_Male/04_Classification/a_script/04_classification.R` |

- **Method:** `arc_classification()` keeps this sex's rows from `09_Classification`: AUC (1 = higher in PAS) with DeLong 95% CI, left empty at perfect separation; exact Wilcoxon p; BH within level; r with PAS beside each row.
- **Report:** perfect separations against all 252 relabellings, AUC distributions, then the ROC curve of every feature at p ≤ 0.05.
- **Result:** Pathway scores: 20 of 889 at p ≤ 0.05 (44 expected); 3 perfect separators (null median 2, p 0.41). Modules 2 of 18, complexes 1 of 93, clusters 0 of 13. None at FDR ≤ 0.05.
- **Note:** Proteins are left out, because their ROC repeats the protein contrast in `01_Proteins`.
