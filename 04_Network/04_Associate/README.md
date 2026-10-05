# 04_Associate

Asks whether each module eigengene tracks each of the six phenotypes within groups, not just across
them.

- **Method:** `lm` per module and phenotype (phenotype z-scored) under three models: pooled slope with group means in the design, females only, and male-minus-female slope difference; BH across the nine modules per phenotype and model. Supplementary: DGCA differential correlation (`DGCA::ddcorAll`, Pearson, BH) for the four comparisons. Pooled module-phenotype correlations (`corPvalueStudent`) on a descriptive sheet.
- **Inputs:** `../01_Build_Modules/c_data/modules.rds`, `02_Differential_Expression/02_Contrasts/c_data/contrasts.rds`, `00_Input/phenotypes.xlsx`
- **Outputs:** `c_data/04_associate.xlsx` (sheets: summary, association, pooled_correlation, dgca, plot_index), `b_reports/04_associate.pdf`
- **Result:** 2 of 162 slopes reach p ≤ 0.05 (chance gives about 8): brown with CSA, pooled (p = 0.014, FDR 0.12), and black with `conductance_sol_pyr`, slope difference (p = 0.048, FDR 0.43). None reaches FDR ≤ 0.05. DGCA: 10 of 216 at p ≤ 0.05, 2 at FDR ≤ 0.05 (pink in PAS_vs_VEH_female and black in sex_in_VEH, both `conductance_sol_pyr`, q 0.038).
- **Note:** the pooled correlations restate group differences (pooling across groups), so they are not tests of association.
