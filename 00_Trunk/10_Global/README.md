# 10_Global

Asks what shapes the whole proteome and the phenotypes: sex, treatment, their interaction, or the pellet's mitochondrial share. Context for the arcs, not one of their questions.

| | |
|---|---|
| Reads | `07_Fit/c_data/fit.rds`, `01_Filtering/c_data/precursors_filtered.rds` |
| Writes | `c_data/10_global.xlsx`, `c_data/global.rds` (read by `Figures/`), `b_reports/10_global.pdf` |
| Run | `Rscript 00_Trunk/10_Global/a_script/10_global.R` |

- **Method:** `prcomp` on the normalised full proteome (centred, not scaled) and on the six phenotypes (scaled); `vegan::adonis2` (Euclidean, 999 permutations, sequential terms) for sex × treatment; `betadisper` by group; `envfit` of the MitoCarta share on PC1-PC2; R² of each factor on PC1 to PC4; each technical variable (precursors, protein groups, raw signal, MitoCarta share) against treatment × sex.
- **Result:** proteome: sex R² 0.38 (p 0.001), treatment 0.04 (p 0.27), interaction 0.04 (p 0.26); dispersion alike (p 0.58); the MitoCarta share fits PC1-PC2 with R² 0.89 (p 0.001). Phenotypes: sex R² 0.18 (p 0.002), treatment 0.11 (p 0.014), interaction 0.07 (p 0.095). Every technical variable differs by sex (p ≤ 0.002) and none by treatment or interaction (p ≥ 0.099).
- **Note:** Sex and pellet composition coincide on PC1, so this step cannot separate them. The arcs compare PAS with vehicle within a sex, where both share that composition.
