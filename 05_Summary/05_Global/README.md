# 05_Global

Asks what shapes the whole proteome: sex, treatment, their interaction, or the mitochondrial
content of the pellet.

- **Method:** `prcomp` on the normalised full universe (20 mice, centred, not scaled); `vegan::adonis2` (Euclidean, 999 permutations, sequential terms) for sex × treatment; `vegan::betadisper` by group with `permutest`; MitoCarta share fitted onto PC1 to PC2 with `vegan::envfit` (999 permutations); R² of sex, treatment, interaction and MitoCarta share on PC1 to PC4.
- **Inputs:** `02_Differential_Expression/02_Contrasts/c_data/contrasts.rds`
- **Outputs:** `c_data/05_global.xlsx` (sheets: scores, variance, permanova, component_r2, mito_vector), `c_data/global.rds` (read by `06_Figures/F01`)
- **Result:** sex explains 38% of the distance (F = 11.1, p = 0.001); treatment 4% (p = 0.27) and the interaction 4% (p = 0.26). Group dispersions are alike (p = 0.58). PC1 (41%) is sex (R² 0.91) and the MitoCarta share (0.86), which coincide; treatment shows only on PC4 (7%, R² 0.23). The MitoCarta vector fits PC1 to PC2 with R² 0.89 (p = 0.001).
- **Note:** sex and MitoCarta share are confounded on PC1, so this stage cannot separate a sex effect from a difference in mitochondrial enrichment.
