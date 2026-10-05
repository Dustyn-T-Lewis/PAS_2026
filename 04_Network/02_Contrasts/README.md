# 02_Contrasts

Tests whether each module eigengene moved with treatment or sex, using the same five contrasts as
the protein analysis.

- **Method:** per module, `lm(eigengene ~ 0 + group)` and the five contrasts by `emmeans::contrast` (df 16); BH across the nine modules within each contrast (too few for limma's empirical Bayes prior). WGCNA module-trait Pearson r with treatment (PAS = 1) and sex (male = 1), `corPvalueStudent`.
- **Inputs:** `../01_Build_Modules/c_data/modules.rds`, `02_Differential_Expression/02_Contrasts/c_data/contrasts.rds`
- **Outputs:** `c_data/02_contrasts.xlsx` (sheets: summary, results, module_trait, plot_index), `c_data/module_contrasts.rds`, `b_reports/02_contrasts_<contrast>.pdf` (5) and `b_reports/02_contrasts_module_trait.pdf`
- **Result:** no module reaches p ≤ 0.05 for PAS_vs_VEH_female, PAS_vs_VEH_male or treatment_by_sex. At p ≤ 0.05 and FDR ≤ 0.05: sex_in_VEH, turquoise (up in males, t = 13.9, FDR 2e-9) and blue (down, t = -9.9); sex_in_PAS, turquoise, blue, brown (FDR 0.009), black (0.016) and pink (0.005). Magenta (in VEH and in PAS) and pink (in VEH) reach p ≤ 0.05 with FDR 0.052 to 0.055.
- **Note:** module-trait r with sex matches: turquoise 0.98, blue -0.96, pink -0.71; every treatment r is |r| ≤ 0.30, all p > 0.19.
