# Evidence: female

Read 2026-10-16 after the coordinator's rerun of run_all.R. `F` = `03_Comparisons/Female/c_data/`.

| Drawn element | Value in the output | Source |
|---|---|---|
| contrast PAS_F − VEH_F | PAS_vs_VEH_female: PAS_F 1, VEH_F −1 | `01_Processing/03_Design/c_data/03_design.xlsx`, sheet `contrasts` |
| read from 02_Statistics | inputs fit.rds, enrichment.rds, classification.rds | `03_Comparisons/code/*.R` input blocks |
| 177 at p ≤ 0.05 (98 expected) | full: down 95 + up 82 = 177; expected 49.90 + 48.20 = 98.1 | `F/01_proteins.xlsx`, sheet `summary`, universe full, `p_le_0.05`, `expected_by_chance`; also `02_Statistics/01_Fit/c_data/01_fit.xlsx` `contrast_summary` (177, 98.10) |
| 3 at FDR ≤ 0.05 | full down 3 + up 0 | `F/01_proteins.xlsx`, `summary`, `fdr_le_0.05` |
| GO:BP 10, Hallmark 1, MitoCarta 2 at FDR ≤ 0.05 | GO:BP 10, Hallmark 1, KEGG 0, Reactome 0, MitoCarta 2 | `F/02_pathways.xlsx`, sheet `summary`, `fdr_le_0.05` |
| muscle contraction up, cell cycle down | FDR ≤ 0.05 GO:BP: muscle system process, (skeletal/striated) muscle contraction NES > 0; cell cycle (NES −1.89, padj 0.014), mitotic cell cycle NES < 0 | `F/02_pathways.xlsx`, sheet `fgsea`, `padj <= 0.05` |
| capillary density +138% | pct_change 138.2 (p 0.0052) | `F/03_phenotypes.xlsx`, sheet `phenotype_effects` |
| CSA +22%, satellite cells +47% | csa_mean 21.6 (p 0.0020); satellite_cells_per_fiber 47.3 (p 0.0093) | same |
| conductance unchanged | conductance_sol_pyr p 0.45, conductance_pla_pyr p 0.30, conductance_sol_oct p 0.35 | same |
| 16 pathway scores (null median 4, p 0.13) | observed 16, null_median 4, p 0.12698 (32/252) | `F/04_classification.xlsx`, sheet `perfect_separations`, level pathway; `02_Statistics/03_Classification/c_data/classification.rds` `$perfect` |
| 1 of 8 modules (null median 0, p 0.063) | features 8, observed 1, null_median 0, p 0.06349 (16/252) | `F/04_classification.xlsx`, sheets `summary` and `perfect_separations`, level module |
| Footer | c_data holds four .xlsx, no .rds | `ls 03_Comparisons/Female/c_data` |
