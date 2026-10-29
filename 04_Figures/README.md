# 04_Figures

Figures that cross the comparisons, one question each. Each main figure folder holds the composite script (`F0N.R`) and one script per panel in `a_script/`, the composite and panel PDFs and PNGs in `b_reports/`, a workbook with one sheet per panel in `c_data/`, and a README with the panel tables. A figure's supplements sit in the same folder: scripts in `a_script/supp/` (`S#.R` and `S#a_….R`), outputs in `b_reports/supp/`, workbooks as `c_data/S#.xlsx`. Every script reads saved outputs only and fits no model.

![Figure 1: design, phenotype effects and structure, proteome structure, protein counts and cross-sex signatures](F01_Global/b_reports/F01.png)

| Figure | Asks | Panels | Supplements |
|---|---|---|---|
| [`F01_Global`](F01_Global/README.md) | What the design is, what PAS did to the phenotypes, and how the proteome is structured | design; phenotype effects and PCA; proteome PCA; proteins per contrast; cross-sex signatures | S1 processing: filtering, run depth, normalisation, surrogate variable, MitoCarta share, composition markers. S2 structure: component R², module heatmap with AUC, perfect separators |
| [`F02_Proteins`](F02_Proteins/README.md) | Which proteins are worth following up in each sex | candidate STRING networks per sex; every mouse for proteins passing three or more routes | S3 protein lists: three Π heatmaps with cluster ORA, volcanoes. S4 candidate routes per sex. S5 sensitivity: counts per variant, baseline proteins across variants |
| [`F03_Pathways`](F03_Pathways/README.md) | Which pathways each comparison moves | one grid: shown sets per contrast (female, male, Treatment × sex) and MitoCarta, drawn in all three contrasts, with distinct-set counts per contrast | S6 pathway sets: GO Slim compartments, literature hypotheses |
| [`F04_Phenotype_pathways`](F04_Phenotype_pathways/README.md) | Which pathways track each phenotype within each sex | one grid: female and male trait blocks, within-treatment model, with distinct-set counts per trait and sex | S7 protein–trait associations: counts at p ≤ 0.05 by direction, top three proteins per trait |
| [`F05_Complexes`](F05_Complexes/README.md) | Which complexes and mitochondrial pathways PAS shifts | CORUM and MitoCarta dot matrix against contrasts and traits | none |

Shared drawing code (theme, palettes, PCA panel, pathway grid and its set-count bars, heatmaps, candidate network, candidate dots, route matrix) lives in `00_Setup/Shared_functions/figures.R`; the barcode is drawn in F01's own panel script.
