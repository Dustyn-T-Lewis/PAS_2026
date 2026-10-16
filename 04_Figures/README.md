# 04_Figures

Figures that cross the comparisons, one question each. Each figure folder holds one script per panel in `a_script/panels/` and a composite script that runs them and stitches the result; panel and composite PDFs and PNGs in `b_reports/`; a workbook with one sheet per panel in `c_data/`; and a README with a panel table. Every script reads saved outputs only and fits no model.

![Figure 1: design, phenotype effects and structure, proteome structure, protein counts and cross-sex signatures](F01_Global/b_reports/F01_composite.png)

| Figure | Asks | Panels |
|---|---|---|
| [`F01_Global`](F01_Global/README.md) | What the design is, what PAS did to the phenotypes, and how the proteome is structured | design; phenotype effects and PCA; proteome PCA; proteins per contrast; cross-sex signatures |
| [`F02_Proteins`](F02_Proteins/README.md) | Which proteins are worth following up in each sex | candidate STRING networks per sex; every mouse for proteins passing three or more routes |
| [`F03_Pathways`](F03_Pathways/README.md) | Which pathways each comparison moves | full-proteome matrix; MitoCarta matrix |
| [`F04_Phenotype_pathways`](F04_Phenotype_pathways/README.md) | Which pathways track each phenotype within each sex | female and male matrices, pooled and within treatment |
| [`F05_Modules`](F05_Modules/README.md) | How the protein modules relate to PAS and to the phenotypes (WGCNA version) | module × trait heatmap with stability; ROC curves; eigengene scatters |
| [`F05_Complexes`](F05_Complexes/README.md) | Which complexes and mitochondrial pathways PAS shifts (alternative Figure 5) | CORUM and MitoCarta dot matrix against contrasts and traits |
| [`S1_Processing`](S1_Processing/README.md) | What processing kept and chose, and what else is in the pellet | filtering; run depth; normalisation; surrogate variable; MitoCarta share; composition markers |
| [`S2_Proteins`](S2_Proteins/README.md) | Each comparison's own protein list | three Π heatmaps; volcanoes |
| [`S3_Pathway_sets`](S3_Pathway_sets/README.md) | The GO Slim compartments and the sets named in advance | compartments; literature hypotheses |
| [`S4_Protein_associations`](S4_Protein_associations/README.md) | How many proteins track each phenotype | counts against chance; strongest protein per trait |
| [`S5_Structure`](S5_Structure/README.md) | Whether perfect separation beats chance; what the PCs track | perfect separators; component R² |
| [`S6_Protein_patterns`](S6_Protein_patterns/README.md) | Which proteins PAS changes in either sex, by response pattern | F-test heatmap; ORA per pattern |
| [`S7_Candidates`](S7_Candidates/README.md) | Which candidate routes each protein passed | route matrix per sex; full ranked table |
| [`S8_Sensitivity`](S8_Sensitivity/README.md) | Whether the protein results survive other processing and model choices | counts per variant; baseline proteins across variants |

Shared drawing code (theme, palettes, pathway matrix, heatmap, barcode, candidate network, route matrix) lives in `00_Setup/Shared_functions/figures.R`.
