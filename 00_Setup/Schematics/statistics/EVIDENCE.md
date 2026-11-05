# Evidence: statistics

Read 2026-10-16 after the coordinator's rerun of run_all.R.

| Drawn element | Value | Source |
|---|---|---|
| Seven steps in number order | 01_Fit ... 07_Sensitivity | `run_all.R` `statistics`; `02_Statistics/README.md` |
| No arrows between cards | not a chain: 02-07 read 01's fit.rds (and 05 reads proteins.rds; 06 reads enrichment and modules; 07 reads fit and processing files) | each step README "Reads" row |
| 01: limpa dpcDE for proteins, limma for scores | | `01_Fit/README.md` Proteins, Scores and eigengenes |
| 01: TREAT and lfsr | `protein_contrasts` columns `treat_p`, `treat_fdr`, `logFC_shrunk`, `lfsr` | `02_Statistics/01_Fit/c_data/01_fit.xlsx` |
| 01: ordinary t for module eigengenes | module rows of `score_contrasts` have df = 15.0 (residual df, no prior df); other levels 18.6 or 19.8 | `02_Statistics/01_Fit/c_data/fit.rds` `$score_contrasts`, `count(level, df)` |
| 02: fgsea on GO:BP, Reactome, Hallmark, KEGG, MitoCarta | collections in sheet `summary` | `02_Statistics/02_Enrichment/c_data/02_enrichment.xlsx`, `summary` |
| 02: camera beside fgsea | columns `camera_p`, `camera_fdr` in sheet `fgsea`; script comment "camera is unweighted" (02_enrichment.R line 202) | `02_enrichment.xlsx` `fgsea`; `02_enrichment.R` |
| 02: fry for cross-sex signatures | sheet `signatures` (fry direction, p, mixed p) | `02_enrichment.xlsx` |
| 03: single-score AUC within each sex; 252 relabellings | `perfect$relabellings` = 252 for every row | `02_Statistics/03_Classification/c_data/classification.rds` `$perfect` |
| 04: sex R² 0.39, treatment R² 0.04 | sex r2 0.3936, treatment r2 0.0380 | `02_Statistics/04_PCA/c_data/04_pca.xlsx`, sheet `permanova` |
| 05: marker-group scores of pellet contents; PAS vs VEH within each sex | groups SR, ER, myofibril, nerve, red cell, vessel, haemoglobin, MitoCarta share, per sex, Welch p | `02_Statistics/05_Composition/c_data/05_composition.xlsx`, sheet `tests` |
| 06: six routes per protein; STRING communities | route columns statistics, every_mouse, both_sexes, phenotype, hub, pathway; column `community` | `02_Statistics/06_Candidates/c_data/06_candidates.xlsx`, sheets `summary`, `candidates` |
| 07: refits under seven changed choices | variants besides baseline: slope 0.5, slope 1.0, run depth instead of sv1, no adjustment, no sample weights, without PAS_M1, pellet share covariate | `02_Statistics/07_Sensitivity/c_data/07_sensitivity.xlsx`, sheet `counts`, `unique(variant)` |
| Five contrasts and context role | PAS_vs_VEH_female, PAS_vs_VEH_male, treatment_by_sex, sex_in_VEH (context), sex_in_PAS (context) | `01_Processing/03_Design/c_data/03_design.xlsx`, sheets `contrasts`, `comparisons` |
| Readers 03_Comparisons and 04_Figures | | `README.md` data-flow paragraph; `run_all.R` |
