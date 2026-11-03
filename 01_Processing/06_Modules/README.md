# 06_Modules

Co-abundance modules by WGCNA, built on proteins detected in every mouse after removing sex and run depth, so the network follows co-variation within each sex.

| | |
|---|---|
| Reads | `01_Processing/02_Quantification/c_data/proteins.rds`, `01_Processing/01_Filtering/c_data/precursors_filtered.rds` (run depth), `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` |
| Writes | `c_data/06_modules.xlsx`, `c_data/modules.rds`, `b_reports/06_modules.pdf` |
| Run | `Rscript 01_Processing/06_Modules/a_script/06_modules.R` |

## Method

- Features: the 827 proteins detected in all 20 mice. limpa gives undetected cells tied values, which would group proteins by when they go missing.
- Adjustment: `limma::removeBatchEffect(batch = sex, covariates = log2 run depth, design = ~treatment)`.
- Network: `pickSoftThreshold` (signed, bicor, R² cut 0.85) and `blockwiseModules` (signed TOM, minModuleSize 20, mergeCutHeight 0.25, deepSplit 2).
- Eigengenes come from the adjusted matrix, so they carry no sex difference; only within-sex contrasts on modules are read.
- Stability: the network is rebuilt on 100 subsamples of 16 mice (four per group); each module gets its median best-match Jaccard.
- Module colours are assigned per build and change if the input changes.

## Results

- Power 20; 8 modules plus 402 unassigned (grey) proteins.
- Turquoise (105 proteins, 84% MitoCarta) is named by oxidative phosphorylation, red (47) by myogenesis.
- Median Jaccard runs from 0.29 (pink) to 0.67 (red).
