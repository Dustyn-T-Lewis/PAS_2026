# 04_PCA

Asks what shapes the whole proteome and the phenotypes: sex, treatment, their interaction or the pellet's mitochondrial share. This is context for the comparisons, not one of their questions.

| | |
|---|---|
| Reads | `02_Statistics/01_Fit/c_data/fit.rds`, `01_Processing/01_Filtering/c_data/precursors_filtered.rds` |
| Writes | `c_data/04_pca.xlsx`, `c_data/pca.rds`, `b_reports/04_pca.pdf` |
| Run | `Rscript 02_Statistics/04_PCA/a_script/04_pca.R` |

## Method

- `prcomp` on the normalised full proteome (centred) and on the six phenotypes (scaled). No surrogate variable.
- `vegan::adonis2` (Euclidean, 999 permutations, sequential terms: sex, treatment, treatment × sex); `betadisper` by group; `envfit` of the MitoCarta share on PC1 and PC2.
- R² of sex, treatment, interaction, MitoCarta share and run depth on PC1 to PC4; each technical variable tested against treatment × sex.

## Results

- Proteome: sex R² 0.39 (p 0.001), treatment 0.04 (p 0.29), interaction 0.04 (p 0.24); dispersion alike (p 0.49). The MitoCarta share fits PC1 and PC2 with R² 0.91.
- Phenotypes: sex R² 0.18 (p 0.002), treatment 0.11 (p 0.014), interaction 0.07 (p 0.095).
- Every technical variable differs by sex (p ≤ 0.002) and none by treatment or interaction (p ≥ 0.099).
- Sex and pellet composition coincide on PC1, so this step cannot separate them.
