# 03 · Pathway Enrichment

Which gene sets each contrast moved, which set scores alone separate two groups, and which sets
track the six phenotypes. The full proteome (1,842 proteins) is tested against Hallmark, Reactome,
KEGG MEDICUS and GO:BP; the mitochondrial universe (568 proteins) against MitoCarta 3.0 pathways;
GO Slim cellular components get their own test in `07_Compartments`. Tests come from
`enrichVolcano::run_enrichment()`: fgsea leads, camera and fry check it. BH runs within each
collection. `dedup_terms(method = "enrichmentmap")` flags redundant terms for display only; no
term is dropped from a test. A term counts as significant at p ≤ 0.05 and FDR ≤ 0.05. Every
workbook opens with a `read_me` sheet and ends with `input_manifest` (files read, md5) and
`package_versions`.

| Sub-stage | Does | Writes |
|---|---|---|
| `00_Gene_Sets` | Loads the collections, GO Slim ancestors, ROS and literature-informed sets | `gene_sets.rds`, `00_gene_sets.xlsx` |
| `01_Scores` | One singscore per mouse per set | `set_scores.rds`, `01_scores.xlsx`, overview PDF |
| `02_Contrasts` | fgsea, camera and fry for all five contrasts | `enrichment.rds`, `02_contrasts.xlsx`, one PDF per contrast |
| `03_Classify` | AUC of each set score for each two-group comparison | `03_classify.xlsx`, PDF |
| `04_Associate` | fgsea on phenotype rankings, limma on singscores | `04_associate.xlsx`, one PDF per phenotype |
| `05_Volcano` | Protein volcano with fgsea terms ringed around it | one PDF per contrast |
| `06_NES_Scatter` | Female against male NES, sex in VEH against sex in PAS | two PDFs |
| `07_Compartments` | fgsea and fry on 25 GO Slim compartments | `compartments.rds`, `compartments_fry.rds`, `07_compartments.xlsx` |
| `08_Display` | Picks the terms the figures draw | `display.rds`, `08_display.xlsx` |

Run order: `00`, `02`, `01`, `05`, `06`, `07`, `08`; then `03` and `04` after
`02_Differential_Expression/04_Associate`.
