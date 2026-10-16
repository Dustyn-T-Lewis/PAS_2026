# 06_Modules

Co-abundance modules by WGCNA's standard workflow, on the proteins detected in every mouse, built after removing sex and measured run depth so the network follows co-variation within each sex.

| | |
|---|---|
| Reads | `01_Processing/02_Quantification/c_data/proteins.rds`, `01_Processing/01_Filtering/c_data/precursors_filtered.rds` (run depth), `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` |
| Writes | `c_data/06_modules.xlsx`, `c_data/modules.rds` (read by 02_Statistics/01_Fit, 02_Statistics/06_Candidates and `04_Figures/F05_Modules`), `b_reports/06_modules.pdf` |
| Run | `Rscript 01_Processing/06_Modules/a_script/06_modules.R` |

- **Features:** the 827 proteins detected in all 20 mice. limpa gives every undetected mouse the same value for a protein, and those ties would group proteins by when they go missing.
- **Adjustment:** `limma::removeBatchEffect(batch = sex, covariates = log2 run depth, design = ~treatment)`. The earlier build removed sex only; its two largest modules tracked run depth (r 0.99 and −0.91 with the surrogate variable within sex).
- **Network:** `pickSoftThreshold` (signed, bicor, R² cut 0.85, powers to 30) and `blockwiseModules` (signed TOM, minModuleSize 20, mergeCutHeight 0.25, deepSplit 2).
- **Eigengenes:** from the same adjusted matrix the modules were built on (Langfelder, Bioconductor support). The adjustment removes sex differences from them, so only within-sex contrasts on modules are read.
- **Stability:** the network is rebuilt on 100 subsamples of 16 mice (four per group); `sizes` gives each module's median best-match Jaccard and quartiles, `stability` every resample.
- **Result:** power 20; 8 modules plus 402 grey proteins. Turquoise (105 proteins, 84% MitoCarta) is named by oxidative phosphorylation, red (47) by myogenesis. Median Jaccard 0.29 (pink) to 0.67 (red).
- **Basis:** Langfelder & Horvath 2008 and the WGCNA FAQ (signed networks, bicor, at least 15 samples per network); Parsana et al. 2019 on removing confounders before building co-expression networks. Consensus or per-sex networks would need more than the 10 mice per sex.
- **Note:** module colours are assigned on each build and change if the input changes.
