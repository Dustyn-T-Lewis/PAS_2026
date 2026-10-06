# 01_Design

Fixes the model and the five contrasts every later stage uses, and checks their signs before
anything is fitted. Nothing it fits is reported.

- **Method:** `model.matrix(~ 0 + group)` (four cell means, no intercept; 20 samples, rank 4, 16 residual df) and `limma::makeContrasts()`. One sample per mouse, so nothing is blocked. An unweighted `lmFit()` checks that each contrast reproduces differences of raw group means.
- **Inputs:** `01_Preprocess/02_Quantification/c_data/proteins.rds`
- **Outputs:** `c_data/01_design.xlsx` (sheets: design, contrasts, contrast_check, roles, groups), `c_data/design.rds` (design, contrasts, primary, supplementary, and `universes`: full and mito logical vectors over 1,842 proteins; read by `02_Contrasts`). No PDF.
- **Result:** Four groups of five mice. Fitted contrasts match hand-computed group-mean differences to 2.8e-14. `treatment_by_sex` equals `sex_in_VEH - sex_in_PAS` to within 1e-12.
- **Note:** The script stops if design rows and matrix columns differ in order, if any contrast fails the group-mean check, or if the interaction sign is wrong. There is no pooled `PAS_vs_VEH` contrast: the female and male PAS effects are uncorrelated across proteins (r = -0.07), so a pooled average would describe neither sex.
