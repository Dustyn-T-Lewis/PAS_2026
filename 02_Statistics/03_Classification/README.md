# 03_Classification

Asks how well each pathway score, module eigengene, complex score and cluster score alone separates PAS from vehicle mice within each sex. Proteins are left out because their ROC repeats the protein contrast.

| | |
|---|---|
| Reads | `02_Statistics/01_Fit/c_data/fit.rds` |
| Writes | `c_data/03_classification.xlsx`, `c_data/classification.rds`, `b_reports/03_classification.pdf` |
| Run | `Rscript 02_Statistics/03_Classification/a_script/03_classification.R` |

## Method

- `pROC::roc` AUC (1 = higher in PAS) with DeLong 95% CI; exact Wilcoxon p; BH within sex and level.
- The count of perfect separators is compared with all 252 relabellings of the same ten mice.
- Scores enter unadjusted, since a single-feature AUC has no model to carry the surrogate variable.

## Results

- Females: 47 of 1,005 pathway scores at p ≤ 0.05; 16 separate perfectly against a null median of 4 (p 0.13). One of 8 modules separates perfectly (null median 0, p 0.063).
- Males: 28 of 1,005 pathway scores at p ≤ 0.05; 5 perfect against a null median of 2 (p 0.33). No module separates perfectly.
- No level beats the relabelling null at p ≤ 0.05, and no feature reaches FDR ≤ 0.05 in either sex.
