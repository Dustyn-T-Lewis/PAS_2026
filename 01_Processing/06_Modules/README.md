# 06_Modules

Co-abundance modules by WGCNA's standard workflow, built after removing sex so the network follows variation within each sex rather than the sex difference that dominates this proteome.

| | |
|---|---|
| Reads | `01_Processing/02_Quantification/c_data/proteins.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` |
| Writes | `c_data/06_modules.xlsx`, `c_data/modules.rds` (read by 02_Statistics/01_Fit), `b_reports/06_modules.pdf` |
| Run | `Rscript 01_Processing/06_Modules/a_script/06_modules.R` |

- **Method:** `limma::removeBatchEffect(batch = sex, design = ~treatment)`, then `pickSoftThreshold` (signed, bicor, R² cut 0.85) and `blockwiseModules` (signed TOM, minModuleSize 20, mergeCutHeight 0.25, deepSplit 2). Eigengenes come from the unadjusted matrix, so `02_Statistics/01_Fit` can test sex and treatment on them. A module is named only by an over-represented Hallmark or Reactome term at FDR ≤ 0.05.
- **Result:** power 12; 18 modules plus 160 grey proteins. Turquoise is largest (332). `sizes` gives each module's MitoCarta share and its r with the pellet's mitochondrial share, overall and within group.
- **Note:** Module colours are assigned on each build and change if the input changes.
