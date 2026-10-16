# F05_Modules

Figure 5 (WGCNA version): how the protein modules relate to PAS and to the phenotypes. The alternative Figure 5, MitoCarta pathways and CORUM complexes, is `F05_Complexes`.

The modules come from the rebuilt network in `01_Processing/06_Modules`: 827 proteins detected in every mouse, with sex and measured run depth removed before building (the earlier build's two largest modules tracked run depth). The result is eight modules plus 402 unassigned proteins. Eigengenes are tested with ordinary least-squares t on 15 residual df rather than eBayes, because eight features cannot estimate a prior.

![Figure 5](b_reports/F05_composite.png)

| | |
|---|---|
| Reads | `02_Statistics/01_Fit/c_data/fit.rds`, `01_Processing/06_Modules/c_data/modules.rds`, `02_Statistics/03_Classification/c_data/classification.rds` |
| Writes | `b_reports/panels/<panel>.pdf/.png`, `b_reports/F05_composite.pdf/.png` (180 × 160 mm), `c_data/F05.xlsx` |
| Run | `Rscript 04_Figures/F05_Modules/a_script/F05_composite.R` (each panel script also runs alone) |

Each panel is its own script in `a_script/panels/`; the composite runs them, adds the letters, stitches them and checks the cell count, the shown ROC modules against `03_classification.xlsx`, and ten mice per scatter. Panels read saved outputs and compute nothing beyond counts, filtering and scaling.

## Panels

| Panel | Shows | Source |
|---|---|---|
| a | Each module eigengene against the three contrasts and each phenotype within each sex, after YvO's module–trait heatmap. Every cell is the r-equivalent of its t, t / √(t² + df), so value and outline come from one test; solid outline FDR ≤ 0.05 over the eight modules, dashed p ≤ 0.05, bold \|r\| ≥ 0.4. No cell reaches FDR ≤ 0.05. Female PAS raises blue at p 0.052; male PAS raises brown (p 0.042); turquoise (OXPHOS) differs between the sexes' PAS effects (p 0.049). Left: module size, MitoCarta share, stability J (median best-match Jaccard over 100 subsamples of four mice per group; pink is least stable at 0.29) and top annotated set. | `fit.rds`: `score_contrasts`, `score_associations`; `modules.rds`: `sizes` |
| b | ROC curves for the modules at exact Wilcoxon p ≤ 0.05 in either sex (five against five). A module lower in PAS is drawn as 1 − its curve and labelled so. | `classification.rds`: `roc`, `module_curves` |
| c | The four strongest module–phenotype pairs (pooled model), mouse by mouse. The label gives Pearson r across the ten points and the slope p, which also adjusts for the surrogate variable; the two can differ. | `fit.rds`: `score_associations`, `scores`, `phenotypes` |

## Workbook

`c_data/F05.xlsx` holds one sheet per panel with the numbers it draws, plus `read_me`, `input_manifest` and `package_versions`.
