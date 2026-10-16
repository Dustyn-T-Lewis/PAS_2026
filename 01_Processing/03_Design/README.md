# 03_Design

Fixes the model and the five contrasts every later step uses. The model has the four group means and one surrogate variable for unwanted variation.

| | |
|---|---|
| Reads | `01_Processing/02_Quantification/c_data/proteins.rds`, `01_Processing/01_Filtering/c_data/precursors_filtered.rds` (run depth) |
| Writes | `c_data/03_design.xlsx`, `c_data/design.rds` (read by 02_Statistics/01_Fit, 02_Enrichment, 06_Candidates and 07_Sensitivity, and `04_Figures/S1_Processing`), `b_reports/03_design.pdf` |
| Run | `Rscript 01_Processing/03_Design/a_script/03_design.R` |

- **Method:** `model.matrix(~ 0 + group)`, four cell means, plus `sv1` from `sva::sva(n.sv = 1)` on the normalised full proteome (full model the four groups, null model an intercept); 15 residual df; `limma::makeContrasts()`. One sample per mouse, so nothing is blocked.
- **Surrogate variable:** the strongest structure left once the group means are removed. It follows run depth (log2 precursors identified): r = 0.85 across mice, 0.86 in females, 0.92 in males. It does not differ by treatment (Welch p 0.46 in females, 0.71 in males), so it removes noise rather than PAS signal. Without it, five male runs that went deeper than the other five split the male proteome: random relabellings of the ten males put about 11% of proteins at p ≤ 0.05, more than the true PAS split did.
- **Contrasts:** `PAS_vs_VEH_female` (PAS_F − VEH_F), `PAS_vs_VEH_male` (PAS_M − VEH_M), `treatment_by_sex` ((PAS_F − VEH_F) − (PAS_M − VEH_M)); the three comparisons in `03_Comparisons` read these. `sex_in_VEH` and `sex_in_PAS` (male − female) are context only.
- **Basis:** limma User's Guide on cell-means designs and contrasts. Leek & Storey 2007 (*PLoS Genet*): unmodelled factors cost power and can add spurious signal "even for well-designed, randomized studies". The sva vignette fixes one surrogate variable for small studies. Čuklina et al. 2021 (*Mol Syst Biol*) on diagnosing technical variation in proteomics before correcting it. Males ran shallower on average, so part of the baseline sex difference is run depth.
