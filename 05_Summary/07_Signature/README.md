# 07_Signature

Asks whether one sex's PAS signature moves the same way in the other sex.

- **Method:** signature = proteins at p ≤ 0.05 for PAS in one sex (full universe), split up and down; each set tested in the other sex's PAS contrast with `limma::fry` (self-contained, precision weights, stored design and contrasts). Share of signature proteins on the same side of the other sex's t ranking.
- **Inputs:** `02_Differential_Expression/02_Contrasts/c_data/contrasts.rds`
- **Outputs:** `c_data/07_signature.xlsx` (sheets: tests, ranks), `c_data/signature.rds` (read by `08_Signature_Phenotype` and `06_Figures/F03`)
- **Result:** neither signature carries over. Female signature (65 up, 87 down) in males: p = 0.54 and 0.64; male signature (60 up, 73 down) in females: p = 0.60 and 0.999. Between 44% and 62% of signature proteins move the same way in the other sex, against 50% by chance.
- **Note:** no carry-over means the responses are independent, not opposite; whether they differ is tested directly by treatment_by_sex (see `04_Sex_Patterns`).
