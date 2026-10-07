# 07_Fit

Every model the arcs read, fitted once on all 20 mice.

| | |
|---|---|
| Reads | `02_Quantification/c_data/proteins.rds`, `03_Design/c_data/design.rds`, `05_Scores/c_data/scores.rds`, `06_Modules/c_data/modules.rds`, `00_Input/phenotypes.xlsx` |
| Writes | `c_data/07_fit.xlsx`, `c_data/fit.rds` (read by 08 to 10 and the arcs), `c_data/study_full/` and `c_data/study_mito/` (enrichVolcano study tables read by 08), `b_reports/07_fit.pdf` |
| Run | `Rscript 00_Trunk/07_Fit/a_script/07_fit.R` |

- **Proteins:** `limpa::dpcDE(sample.weights = TRUE)` per universe, `contrasts.fit()`, `eBayes()`; BH within universe and contrast; Π = p^|log2 FC| beside p and FDR; r with group (first group = 1) beside each two-group contrast.
- **Scores and eigengenes:** plain limma on the same design and contrasts.
- **Phenotype associations:** each feature regressed on one z-scored phenotype across the ten mice of one sex, PAS and vehicle together (`~0 + sex + z:female + z:male`). This is co-movement, not a PAS effect; r pooled and within each group sits beside it. A refit with the pellet's mitochondrial share (centred within group) is kept as a check.
- **Phenotype effects:** PAS against the same sex's vehicle as % change (Welch t on log values); the interaction from `lm(log(value) ~ treatment * sex)`, as the female response relative to the male.
- **Result (full proteome, 1,842 tests, 92 expected at p ≤ 0.05):** female 152 at p ≤ 0.05, 109 at Π ≤ 0.05, 4 at FDR ≤ 0.05; male 133, 78, 0; interaction 156, 163, 0. Phenotype effects in females: capillary density +138% (p 0.005), CSA +22% (p 0.002), satellite cells +47% (p 0.009); none in males; interaction only for capillary density (p 0.014).
- **Note:** Protein contrasts are identical to v1. The script stops unless the joint association model equals a regression within one sex.
