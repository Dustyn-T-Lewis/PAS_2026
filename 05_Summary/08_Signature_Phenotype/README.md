# 08_Signature_Phenotype

Asks whether PAS changes muscle phenotype in each sex, and whether each sex's PAS signature score
tracks phenotype within groups.

- **Method:** (1) PAS minus VEH per phenotype and sex, in SD units of all 20 mice, Welch `t.test` with 95% CI; interaction p from `lm(z ~ treatment * sex)`. (2) Signature score per mouse, `singscore::simpleScore` on the up and down sets from `07_Signature`. (3) `lm` of score on phenotype under three models: pooled slope with group means in the design, females only, male-minus-female slope difference (36 slopes).
- **Inputs:** `02_Differential_Expression/02_Contrasts/c_data/contrasts.rds`, `../07_Signature/c_data/signature.rds`, `00_Input/phenotypes.xlsx`
- **Outputs:** `c_data/08_signature_phenotype.xlsx` (sheets: phenotype_effects, scores, slopes), `c_data/signature_phenotype.rds` (read by `06_Figures/F04`)
- **Result:** in females PAS raises capillary density (+1.7 SD, p = 0.019; interaction p = 0.014), CSA (+0.7 SD, p = 0.003) and satellite cells per fibre (+1.3 SD, p = 0.015). No male phenotype and no conductance changes at p ≤ 0.05. Two of 36 slopes reach p ≤ 0.05 (about 2 by chance), both the male signature against `conductance_sol_oct` (females only p = 0.001; slope difference p = 0.018).
- **Note:** the signatures were chosen on these mice, so the slopes are exploratory. Phenotypes come from soleus, plantaris and IHC (muscle to confirm); the proteome is gastrocnemius.
