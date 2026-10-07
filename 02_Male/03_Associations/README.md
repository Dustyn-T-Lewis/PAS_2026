# 03_Associations

Two questions kept apart. Phenotype effects: what PAS did to each phenotype in males (a treatment effect). Phenotype associations: which proteins and scores move with each phenotype across the males' ten PAS and vehicle mice (co-movement, not a PAS effect).

| | |
|---|---|
| Reads | `00_Trunk/07_Fit/c_data/fit.rds`, `00_Trunk/08_Enrichment/c_data/enrichment.rds` |
| Writes | `c_data/03_associations.xlsx`, `c_data/associations.rds`, `b_reports/03_associations.pdf` |
| Run | `Rscript 02_Male/03_Associations/a_script/03_associations.R` |

- **Method:** `arc_associations()` keeps this sex's rows from `07_Fit`: % change with Welch t on log values; each feature regressed on one z-scored phenotype across the ten mice (limma), with r pooled, in vehicle and in PAS beside it; a refit with the pellet's mitochondrial share as the one check. fgsea on the protein association t comes from `08_Enrichment`.
- **Report:** phenotype effects, association counts by level and phenotype, the mito-share check, r heatmaps of modules, complexes and clusters against treatment and each phenotype, the strongest robust protein per phenotype (keeps sign and p ≤ 0.05 with any one mouse left out), then a scatter grid of every association at p ≤ 0.05.
- **Result:** Phenotype effects: none at p ≤ 0.05; the closest is satellite cells per fibre, +43% (p 0.086). Phenotype associations at p ≤ 0.05: 2,185 of 11,052 protein-phenotype tests in the full proteome (553 expected), 488 at FDR ≤ 0.05; full-proteome pathway scores 167 of 5,076, fewer than the 254 expected; complexes 21 of 558; modules 4 of 108; clusters 5 of 78.
- **Note:** Each association grid page states the phenotype effect of PAS in this sex; where it is p > 0.05, the associations on that page follow variation between mice, not a PAS response. The fgsea on association t permutes genes; no set's singscore association reaches FDR ≤ 0.05.
