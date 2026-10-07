# 03_Design

Fixes the model and the five contrasts every later step uses, and checks their signs before anything is fitted.

| | |
|---|---|
| Reads | `02_Quantification/c_data/proteins.rds` |
| Writes | `c_data/03_design.xlsx`, `c_data/design.rds` (read by 07). No PDF. |
| Run | `Rscript 00_Trunk/03_Design/a_script/03_design.R` |

- **Method:** `model.matrix(~ 0 + group)`, four cell means, 16 residual df; `limma::makeContrasts()`. One sample per mouse, so nothing is blocked; no covariates (run order was randomised).
- **Contrasts:** `PAS_vs_VEH_female` (PAS_F − VEH_F), `PAS_vs_VEH_male` (PAS_M − VEH_M), `treatment_by_sex` ((PAS_F − VEH_F) − (PAS_M − VEH_M)); the arcs read these three. `sex_in_VEH` and `sex_in_PAS` (male − female) are context only.
- **Note:** The script stops if design rows and matrix columns differ in order, if a contrast does not reproduce raw group-mean differences, or if the interaction sign is wrong.
