# Evidence: sex_difference

Read 2026-10-16 after the rerun of run_all.R; every number re-checked 2026-10-17 against the c_data workbooks at HEAD 4a15743 (unchanged). `S` = `03_Comparisons/Sex_Difference/c_data/`.

| Drawn element | Value in the output | Source |
|---|---|---|
| contrast (PAS_F − VEH_F) − (PAS_M − VEH_M) | treatment_by_sex: PAS_F 1, PAS_M −1, VEH_F −1, VEH_M 1 | `01_Processing/03_Design/c_data/03_design.xlsx`, sheet `contrasts` |
| 162 at p ≤ 0.05 | full: down 59 + up 103 = 162 | `S/01_proteins.xlsx`, sheet `summary`, universe full, `p_le_0.05` (re-read 2026-10-17 at HEAD 4a15743) |
| none at FDR ≤ 0.05 | 0 | `S/01_proteins.xlsx`, `summary`, `fdr_le_0.05` |
| GO:BP 20, Reactome 6, Hallmark 3, KEGG 2, MitoCarta 11 at FDR ≤ 0.05 | as drawn | `S/02_pathways.xlsx`, sheet `summary`, `fdr_le_0.05` |
| muscle contraction and extracellular matrix higher in the female response | NES > 0 at padj ≤ 0.05: GOBP muscle contraction, REACTOME striated muscle contraction, muscle contraction, extracellular matrix organization, degradation of the extracellular matrix | `S/02_pathways.xlsx`, sheet `fgsea`; positive = PAS raised it more in females (`Sex_Difference/README.md` first paragraph, matching the contrast signs) |
| phenotype interaction: capillary density only (p 0.014) | capillary_density P.Value 0.0136; the other five 0.22 to 0.95 | `02_Statistics/01_Fit/c_data/01_fit.xlsx`, sheet `phenotype_interaction` |
| no classification step | only 01_proteins.xlsx and 02_pathways.xlsx | `ls 03_Comparisons/Sex_Difference/c_data`; `03_Comparisons/README.md` table |
| Footer | c_data holds two .xlsx, no .rds | `ls 03_Comparisons/Sex_Difference/c_data` |
