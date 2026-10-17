# S4_Protein_associations

Supplementary Figure 4: how many proteins track each phenotype, and the strongest one.

![Supplementary Figure 4](b_reports/S4_composite.png)

| | |
|---|---|
| Reads | `02_Statistics/01_Fit/c_data/fit.rds` |
| Writes | `b_reports/panels/<panel>.pdf/.png`, `b_reports/S4_composite.pdf/.png` (180 × 135 mm), `c_data/S4.xlsx` |
| Run | `Rscript 04_Figures/S4_Protein_associations/a_script/S4_composite.R` (each panel script also runs alone) |

Each panel is its own script in `a_script/panels/`; the composite runs them, adds the letters, stitches them and writes one workbook. Panels read saved outputs and compute nothing beyond counts, filtering and scaling.

## Panels

### a. Association counts

![Panel a](b_reports/panels/a_counts.png)

Each bar counts the full-proteome proteins whose association with one phenotype, within one sex, reaches p ≤ 0.05 under one model (`fit.rds`: `association_summary`): pooled across PAS and vehicle mice (dark grey), pooled with the pellet's MitoCarta share added (light grey), and within treatment (red). The within-treatment model removes the group means and asks whether a protein tracks the trait among mice given the same treatment.

Female satellite cells per fibre stand out under the pooled models (231 and 241), and fall to 156 within treatment. Female capillary density rises from 110 pooled to 169 within treatment, and female soleus pyruvate rises from 67 pooled to 124 within treatment. In males, plantaris pyruvate gives the highest counts under all three models (184, 187 and 164), and satellite cells drop to 47 within treatment. FDR ≤ 0.05 is rare, at most 8 proteins in any bar (male soleus pyruvate, pooled + MitoCarta share). Proteins co-vary and each sex has ten mice, so a count is not itself a test.

### b. Strongest protein per phenotype

![Panel b](b_reports/panels/b_top_proteins.png)

Each row is the protein with the smallest within-treatment p for one phenotype in one sex (`fit.rds`: `protein_associations`). Its two points are the r-equivalent of the limma t, t / sqrt(t² + df), under the pooled (grey) and within-treatment (red) models, joined by a line. Filled points pass FDR ≤ 0.05 and open points do not; where the two estimates coincide, the red point hides the grey one.

Within-treatment r-equivalents run from 0.72 to 0.93 in size, but only three reach FDR ≤ 0.05: Cpa3 with female soleus octanoyl conductance (−0.93, FDR 0.012), Fabp4 with male plantaris pyruvate (0.84, 0.024) and Ufl1 with male soleus pyruvate (−0.83, 0.039). All three also pass pooled. Some rows differ sharply between models: Fhl3 with female capillary density is −0.78 within treatment and 0.02 pooled, and Hnrnpa2b1 with male satellite cells is 0.72 against 0.12. Ssbp1 is the top male protein for both capillary density and CSA. Each row was picked as the best of 1,962, so its within-treatment estimate runs optimistic.

The models come from `02_Statistics/01_Fit`: both regress each protein on one z-scored phenotype within one sex with the surrogate variable fitted within sex; the MitoCarta share in the check model is centred within group. Both models measure co-movement, not a PAS effect.

## Workbook

`c_data/S4.xlsx` holds one sheet per panel with the numbers it draws, plus `read_me`, `input_manifest` and `package_versions`.
