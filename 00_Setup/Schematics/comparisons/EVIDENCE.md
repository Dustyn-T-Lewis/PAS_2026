# Evidence: comparisons

Read 2026-10-16 after the coordinator's rerun of run_all.R.

| Drawn element | Value | Source |
|---|---|---|
| Four steps written once in code/ | 01_proteins.R, 02_pathways.R, 03_phenotypes.R, 04_classification.R | `ls 03_Comparisons/code` |
| Inputs 02_Statistics (fit, enrichment, classification) and 01_Processing (gene sets) | 01 reads fit.rds, enrichment.rds; 02 reads enrichment.rds, gene_sets.rds; 03 reads fit.rds, enrichment.rds; 04 reads classification.rds, fit.rds | `inputs <- c(...)` blocks in each `03_Comparisons/code/*.R` |
| 01 proteins: contrast, tiers by FDR, Π and p | sheet `results` column tier; `summary` columns pi_le_0.05, fdr_le_0.05, p_le_0.05 | `03_Comparisons/Female/c_data/01_proteins.xlsx`, `read_me` and `summary` |
| 02 pathways: fgsea beside camera and singscore | `fgsea` sheet columns NES, p, padj, score_logFC, score_p, camera_p, camera_fdr | `03_Comparisons/Female/c_data/02_pathways.xlsx`, sheet `fgsea` |
| 03 phenotypes: effects and associations | sheets `phenotype_effects`, `protein_associations`, `score_associations` | `03_Comparisons/Female/c_data/03_phenotypes.xlsx` |
| 04 classification: single-score AUC, 252 relabellings | `perfect_separations$relabellings` = 252 | `03_Comparisons/Female/c_data/04_classification.xlsx`, sheet `perfect_separations` |
| Contrasts | PAS_F − VEH_F; PAS_M − VEH_M; (PAS_F − VEH_F) − (PAS_M − VEH_M) | `01_Processing/03_Design/c_data/03_design.xlsx`, sheet `contrasts` |
| Steps run per folder | Female and Male: 01-04 workbooks; Sex_Difference: 01_proteins.xlsx, 02_pathways.xlsx only | `ls 03_Comparisons/*/c_data` |
| Footer (no .rds) | each c_data holds only <step>.xlsx; no code/ script calls saveRDS | `ls 03_Comparisons/*/c_data`; `grep saveRDS 03_Comparisons/code/*.R` (no hit) |
