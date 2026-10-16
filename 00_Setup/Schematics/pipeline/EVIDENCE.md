# Evidence: pipeline

Repo state read 2026-10-16, branch refactor/hrvlr-layout, after the coordinator's full rerun of run_all.R.

| Drawn element | Source |
|---|---|
| Run order 01_Processing -> 02_Statistics -> 03_Comparisons -> 04_Figures; chip "Rscript run_all.R" | `run_all.R`, vectors `processing`, `statistics`, `comparisons`, `figures` |
| 00_Setup inputs: report.parquet, report_stats.xlsx, metadata.csv, phenotypes.xlsx, MitoCarta 3.0 (`Mouse.MitoCarta3.0.xls`), CORUM 5.3 (`corum_complete.txt`), STRING v12 (`STRING_v12/`) | `ls 00_Setup/Input`; `00_Setup/Input/README.md` table. `VEH_PAS_HLI_Master Data.xlsx` is read by nothing (README), so not drawn |
| Shared_functions | `00_Setup/Shared_functions/` (helpers.R, figures.R, comparisons.R) |
| Six processing steps | `run_all.R` `processing`; folder names `01_Processing/0N_*` (04 Gene sets = `04_Gene_Sets`) |
| 1,962 proteins, 595 in MitoCarta | `01_Processing/02_Quantification/c_data/proteins.rds`: `dim(proteins)` = 1962 x 20; `dim(E_norm$mito)` = 595 x 20 |
| Seven statistics steps | `run_all.R` `statistics`; folders `02_Statistics/01_Fit` ... `07_Sensitivity` |
| "run once on all 20 mice" | `02_Statistics/README.md` first line; 20 columns in `proteins.rds` |
| Three comparisons and contrasts | `run_all.R` `comparisons`; `01_Processing/03_Design/c_data/03_design.xlsx` sheets `contrasts` (PAS_vs_VEH_female = PAS_F − VEH_F; PAS_vs_VEH_male = PAS_M − VEH_M; treatment_by_sex = PAS_F − PAS_M − VEH_F + VEH_M) and `comparisons` |
| "four steps written once in code/" | `ls 03_Comparisons/code`: 01_proteins.R, 02_pathways.R, 03_phenotypes.R, 04_classification.R |
| 04_Figures: F01_Global, F02_Proteins, F03_Pathways, F04_Phenotype_pathways, F05_Modules, F05_Complexes, S1–S8 | `run_all.R` `figures` (14 composites; S1_Processing to S8_Sensitivity) |
| Arrows into 03_Comparisons and 04_Figures from 02_Statistics; subtitle "drawn from saved results only" | `README.md` "Data flows one way"; `04_Figures/README.md` ("Every script reads saved outputs only"). Figures also read some 01_Processing files (F02, F05_Modules, F05_Complexes, S1 read proteins.rds, precursors_filtered.rds, design.rds, gene_sets.rds or modules.rds, by grep of `04_Figures/*/a_script`), so the subtitle does not say "statistics only" |
| Footer: a_script/, b_reports/, c_data/, README.md | `README.md` section "Every step" |
