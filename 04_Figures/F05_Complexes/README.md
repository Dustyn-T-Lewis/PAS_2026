# F05_Complexes

Figure 5 (complexes version): CORUM complexes and MitoCarta pathways that PAS shifts, against the contrasts and the phenotypes. The WGCNA version is `F05_Modules`; one of the two goes forward.

![Figure 5, complexes version](b_reports/F05_complexes_composite.png)

| | |
|---|---|
| Reads | `02_Statistics/01_Fit/c_data/fit.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` |
| Writes | `b_reports/panels/a_set_matrix.pdf/.png`, `b_reports/F05_complexes_composite.pdf/.png` (180 × 120 mm), `c_data/F05_complexes.xlsx` |
| Run | `Rscript 04_Figures/F05_Complexes/a_script/F05_complexes_composite.R` (the panel script also runs alone) |

Each set has one singscore per mouse (`01_Processing/05_Scores`), the mean within-mouse rank of its members, so it uses no group label. `01_Fit` tests the scores with limma: the three contrasts, and each trait's slope within each sex (pooled model, adjusted for the surrogate variable). Unlike fgsea in F03, the score test asks whether the set's members move together across mice, not whether they sit high in a protein ranking.

## Panels

### a. Complexes and mitochondrial pathways

![Panel a](b_reports/panels/a_set_matrix.png)

One dot per set and column, for the 21 CORUM complexes and 5 MitoCarta pathways whose singscore differs at p ≤ 0.05 in at least one contrast (`score_contrasts`, `score_associations`, `membership`). The first three columns are the contrasts (PAS against vehicle in females and in males, and Treatment × sex, the female PAS effect minus the male); the rest are each trait's slope within one sex. Fill is the r-equivalent of the limma t, brown when the score is lower in PAS or falls with the trait and teal when it is higher or rises; size is −log10 p, capped at 4; a thick ring marks BH FDR ≤ 0.05 over every set scored in that universe (1,066 full, 44 MitoCarta); the number in brackets counts the set's measured members.

The DDB2 complex is the only set to pass FDR: lower in female PAS (r −0.89, FDR 9.4e-05), different between the sexes' PAS effects (r −0.81, FDR 0.0096), and lower with more female capillaries (r −0.84, FDR 0.0066). Its six measured members are Cops3, Cops4, Cops8, Gps1, Ddb1 and Rbx1, so the measured set is mostly the COP9 signalosome, whose Cops3 and Cops8 also pass FDR ≤ 0.05 at protein level and stay at FDR ≤ 0.10 under every sensitivity variant (S8). At p ≤ 0.05, 21 sets move in female PAS: 13 lower, among them the proteasome, chaperonin, TIM23 and MITRAC complexes, and 8 higher, the oligosaccharyltransferase complexes and the MitoCarta ROS, detoxification and protein homeostasis pathways among them. Only 4 move in male PAS, all lower (two dystrophin–sarcoglycan complexes, Zfand2b, Fe-S cluster biosynthesis). Of the 156 trait cells in each sex, 45 reach p ≤ 0.05 in females and 2 in males.

Complexes share subunits (the 20S core subunits recur across the proteasome entries), so neighbouring rows are not independent results. Rows were chosen on their uncorrected contrast p, so read the female block as a consistent direction across related complexes rather than as 21 separate findings.

## Workbook

`c_data/F05_complexes.xlsx` holds one sheet per panel with the numbers it draws, plus `read_me`, `input_manifest` and `package_versions`.
