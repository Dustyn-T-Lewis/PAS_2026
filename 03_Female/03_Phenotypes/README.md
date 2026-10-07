# 03_Phenotypes

Two questions kept apart. Phenotype effects: what PAS did to each phenotype in females (a treatment effect). Phenotype associations: which proteins and scores move with each phenotype across the females' ten PAS and vehicle mice (co-movement, not a PAS effect).

| | |
|---|---|
| Reads | `02_Statistics/01_Fit/c_data/fit.rds`, `02_Statistics/02_Enrichment/c_data/enrichment.rds` |
| Writes | `c_data/03_phenotypes.xlsx`, `c_data/phenotypes.rds`, `b_reports/03_phenotypes.pdf` |
| Run | `Rscript 03_Female/03_Phenotypes/a_script/03_phenotypes.R` |

- **Method:** `report_phenotypes()` keeps this sex's rows from `02_Statistics/01_Fit`: % change with Welch t on log values; each feature regressed on one z-scored phenotype across the ten mice (limma), with r pooled, in vehicle and in PAS beside it; a refit with the pellet's mitochondrial share as the one check. fgsea on the protein association t comes from `02_Statistics/02_Enrichment`.
- **Report:** phenotype effects, association counts by level and phenotype, the mito-share check, r heatmaps of modules, complexes and clusters against treatment and each phenotype, the strongest robust protein per phenotype (keeps sign and p ≤ 0.05 with any one mouse left out), then a scatter grid of every association at p ≤ 0.05.
- **Result:** Phenotype effects: capillary density +138% (p 0.005), CSA +22% (p 0.002), satellite cells per fibre +47% (p 0.009); the three respiratory conductances do not change (p ≥ 0.30). Phenotype associations at p ≤ 0.05: 803 of 11,052 protein-phenotype tests in the full proteome (553 expected), 20 at FDR ≤ 0.05, 14 of them with satellite cells; full-proteome pathway scores 248 of 5,076 (254 expected), none at FDR; complexes 44 of 558, one at FDR; modules 7 of 108; clusters 5 of 78.
- **Note:** Each association grid page states the phenotype effect of PAS in this sex; where it is p > 0.05, the associations on that page follow variation between mice, not a PAS response. The fgsea on association t permutes genes; no set's singscore association reaches FDR ≤ 0.05.
