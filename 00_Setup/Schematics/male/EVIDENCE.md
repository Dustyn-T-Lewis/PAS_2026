# Evidence: male

Read 2026-10-16 after the coordinator's rerun of run_all.R. `M` = `03_Comparisons/Male/c_data/`.

| Drawn element | Value in the output | Source |
|---|---|---|
| contrast PAS_M − VEH_M | PAS_vs_VEH_male: PAS_M 1, VEH_M −1 | `01_Processing/03_Design/c_data/03_design.xlsx`, sheet `contrasts` |
| 193 at p ≤ 0.05 (98 expected) | full: down 119 + up 74 = 193; expected 54.75 + 43.35 = 98.1 | `M/01_proteins.xlsx`, sheet `summary`; `01_fit.xlsx` `contrast_summary` (193, 98.10) |
| none at FDR ≤ 0.05 | full 0 + 0 | `M/01_proteins.xlsx`, `summary`, `fdr_le_0.05` |
| GO:BP 35, Hallmark 1, KEGG 5 at FDR ≤ 0.05 | GO:BP 35, Hallmark 1, KEGG 5, Reactome 0, MitoCarta 0 | `M/02_pathways.xlsx`, sheet `summary`, `fdr_le_0.05` |
| actin and cytoskeleton organisation down | 34 of 35 GO:BP at FDR ≤ 0.05 have NES < 0, led by supramolecular fiber organization, cytoskeleton organization, actin filament organization | `M/02_pathways.xlsx`, sheet `fgsea`, `padj <= 0.05` |
| no phenotype effect at p ≤ 0.05 | smallest P.Value 0.086 (satellite cells) | `M/03_phenotypes.xlsx`, sheet `phenotype_effects` |
| 5 pathway scores (null median 2, p 0.33) | observed 5, null_median 2, p 0.3333 (84/252) | `M/04_classification.xlsx`, sheet `perfect_separations`, level pathway; `classification.rds` `$perfect` |
| 0 of 8 modules | module features 8, observed 0 | `M/04_classification.xlsx`, sheets `summary` and `perfect_separations` |
| Footer | c_data holds four .xlsx, no .rds | `ls 03_Comparisons/Male/c_data` |
