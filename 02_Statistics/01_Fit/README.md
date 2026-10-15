# 01_Fit

Every model the comparisons read, fitted once on all 20 mice.

| | |
|---|---|
| Reads | `01_Processing/02_Quantification/c_data/proteins.rds`, `01_Processing/03_Design/c_data/design.rds`, `01_Processing/05_Scores/c_data/scores.rds`, `01_Processing/06_Modules/c_data/modules.rds`, `00_Setup/Input/phenotypes.xlsx` |
| Writes | `c_data/01_fit.xlsx`, `c_data/fit.rds` (read by 02_Enrichment to 04_PCA, `03_Comparisons` and `04_Figures/`), `c_data/study_full/` and `c_data/study_mito/` (enrichVolcano study tables read by 02_Enrichment), `b_reports/01_fit.pdf` |
| Run | `Rscript 02_Statistics/01_Fit/a_script/01_fit.R` |

- **Proteins:** `limpa::dpcDE(sample.weights = TRUE)` per universe, `contrasts.fit()`, `eBayes()`; BH within universe and contrast; Π = p^|log2 FC| beside p and FDR; r with group (first group = 1) beside each two-group contrast.
- **Scores and eigengenes:** plain limma on the same design and contrasts.
- **Phenotype associations:** each feature regressed on one z-scored phenotype across the ten mice of one sex, PAS and vehicle together (`~0 + sex + z:female + z:male`). This is co-movement, not a PAS effect; r pooled and within each group sits beside it. A refit with the pellet's mitochondrial share (centred within group) is kept as a check.
- **Phenotype effects:** PAS against the same sex's vehicle as % change (Welch t on log values); the interaction from `lm(log(value) ~ treatment * sex)`, as the female response relative to the male.
- **Result (full proteome, 1,842 tests, 92 expected at p ≤ 0.05):** female 152 at p ≤ 0.05, 109 at Π ≤ 0.05, 4 at FDR ≤ 0.05; male 133, 78, 0; interaction 156, 163, 0. Phenotype effects in females: capillary density +138% (p 0.005), CSA +22% (p 0.002), satellite cells +47% (p 0.009); none in males; interaction only for capillary density (p 0.014).
- **Workbook:** `protein_contrasts` and `score_contrasts` hold every contrast; `protein_associations` and `score_associations` every phenotype association; `association_summary` the counts, with the mito-share check as `adjusted = TRUE`.
- **Note:** The script stops unless the joint association model equals a regression within one sex, and unless every r joins exactly once.
