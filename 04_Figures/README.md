# 04_Figures

The main figures and their supplements, drawn from saved results only. Each folder holds one script per panel plus the composite script (`a_script/`, supplements in `a_script/supp/`), the PDFs and PNGs (`b_reports/`) and one workbook per figure with a sheet per panel (`c_data/`). Drawing code used by more than one figure is in `00_Setup/Shared_functions/figures.R`; code shared by one figure's panels sits beside them in its `a_script/` (`pca_panel.R` in F01, `candidates.R` and `supp/s3_heatmap.R` in F02). A figure with one panel saves only the composite.

| Folder | Holds |
|---|---|
| [`F01_Global`](F01_Global/README.md) | Figure 1 (design, phenotype effects, global structure), S1 (processing), S2 (structure) |
| [`F02_Proteins`](F02_Proteins/README.md) | Figure 2 (candidate proteins per sex), S3 (protein heatmaps and volcanoes), S4 (sensitivity) |
| [`F03_Pathways`](F03_Pathways/README.md) | Figure 3 (pathways per contrast), S5 (compartments and literature sets) |
| [`F04_Phenotype_pathways`](F04_Phenotype_pathways/README.md) | Figure 4 (pathways per phenotype), S6 (protein associations) |
| [`F05_Complexes`](F05_Complexes/README.md) | Figure 5 (complex and MitoCarta scores), no supplement |
