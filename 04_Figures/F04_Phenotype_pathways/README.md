# F04_Phenotype_pathways

Figure 4: within each sex, which pathways track each phenotype across the ten mice.

![Figure 4](b_reports/F04_composite.png)

| | |
|---|---|
| Reads | `02_Statistics/02_Enrichment/c_data/enrichment.rds` |
| Writes | `b_reports/panels/<panel>.pdf/.png`, `b_reports/F04_composite.pdf/.png` (180 × 150 mm), `c_data/F04.xlsx` |
| Run | `Rscript 04_Figures/F04_Phenotype_pathways/a_script/F04_composite.R` (each panel script also runs alone) |

Each panel is its own script in `a_script/panels/`; the composite runs them, adds the letters, stitches them and checks that every row is a kept set at FDR ≤ 0.05 in that sex, listed once, at most four per trait, and drawn in every column. Panels read saved outputs and compute nothing beyond counts, filtering and scaling.

## Panels

| Panel | Shows | Source |
|---|---|---|
| a | Females. fgsea on each protein's association t with each phenotype, pooled over treatment and within treatment. Phenotype by phenotype, the four kept sets at FDR ≤ 0.05 under either model with the smallest FDR that no earlier phenotype already lists; each set drawn for every phenotype under both models. | `enrichment.rds`: `association_fgsea` |
| b | Males, selected and drawn the same way. | `enrichment.rds`: `association_fgsea` |

These are co-movement across mice, not PAS effects. The gene-permutation FDR runs optimistic (Wu & Smyth 2012). camera was added for the contrasts (F03) but not for these associations; the check here is each set's singscore slope on the trait, tested across mice (`score_p`, `score_fdr` in the `association_fgsea` sheet).

## Workbook

`c_data/F04.xlsx` holds one sheet per panel with the numbers it draws, plus `read_me`, `input_manifest` and `package_versions`.
