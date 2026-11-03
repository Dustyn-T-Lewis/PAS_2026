# 03_Design

Fixes the model and the five contrasts every later step uses: four group means plus one surrogate variable for unwanted variation.

| | |
|---|---|
| Reads | `01_Processing/02_Quantification/c_data/proteins.rds`, `01_Processing/01_Filtering/c_data/precursors_filtered.rds` (run depth) |
| Writes | `c_data/03_design.xlsx`, `c_data/design.rds`, `b_reports/03_design.pdf` |
| Run | `Rscript 01_Processing/03_Design/a_script/03_design.R` |

## Method

- `model.matrix(~ 0 + group)` plus `sv1` from `sva::sva(n.sv = 1)` on the normalised full proteome; 15 residual df. One sample per mouse, so nothing is blocked.
- Contrasts: `PAS_vs_VEH_female` (PAS_F − VEH_F), `PAS_vs_VEH_male` (PAS_M − VEH_M) and `treatment_by_sex` ((PAS_F − VEH_F) − (PAS_M − VEH_M)) are the three questions. `sex_in_VEH` and `sex_in_PAS` (male − female) are context only.

## Results

- `sv1` follows run depth (log2 precursors identified): r = 0.85 across mice, 0.86 in females, 0.92 in males.
- It does not differ by treatment (Welch p 0.46 in females, 0.71 in males), so it removes run-depth noise and not PAS signal.
