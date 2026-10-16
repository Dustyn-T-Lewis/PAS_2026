# Evidence: processing

Read 2026-10-16 after the coordinator's rerun of run_all.R.

| Drawn element | Value | Source |
|---|---|---|
| Six steps, in number order | 01 Filtering ... 06 Modules | `run_all.R` `processing`; `01_Processing/README.md` table |
| No arrows between cards | steps do not form a chain: 04 reads 02, not 03; 05 reads 02 and 04; 06 reads 01, 02 and 04 | each step README "Reads" row |
| Input card "00_Setup/Input" | report.parquet, metadata.csv, report_stats.xlsx, MitoCarta, CORUM, STRING | step READMEs, "Reads" rows |
| 01: q ≤ 0.01; contaminants and shared peptides removed | | `01_Processing/01_Filtering/README.md` Method |
| 01: 21,600 precursors, 2,094 protein groups | precursors 21600, protein_groups 2094 | `01_Processing/01_Filtering/c_data/01_filtering.xlsx`, sheet `filter_log`, row "contaminants removed" |
| 02: limpa detection-probability curve | dpcQuant | `02_Quantification/README.md` Method; `02_quantification.xlsx` sheet `dpc_parameters` |
| 02: 1,962 proteins, 595 in MitoCarta | 1962 x 20; mito 595 x 20 | `01_Processing/02_Quantification/c_data/proteins.rds`: `dim($proteins)`, `dim($E_norm$mito)` |
| 03: four group means and one surrogate variable; five contrasts | columns PAS_F, PAS_M, VEH_F, VEH_M, sv1; 5 contrasts | `01_Processing/03_Design/c_data/03_design.xlsx`, sheets `design`, `contrasts` |
| 04: GO:BP, Reactome, Hallmark, KEGG LEGACY, MitoCarta; GO Slim compartments, CORUM, STRING clusters | collections with tested sets: Compartments 25, CORUM 98, GO:BP 755, Hallmark 34, KEGG 33, MitoCarta 44, Reactome 140, STRING 7 | `01_Processing/04_Gene_Sets/c_data/04_gene_sets.xlsx`, sheet `catalog`, `table(collection, tested)`; KEGG LEGACY wording from `04_Gene_Sets/README.md` |
| 05: 1,110 scores per mouse | 1110 rows (full 1066, mito 44): CORUM 98, GO:BP 755, Hallmark 34, KEGG 33, MitoCarta 44, Reactome 139, STRING 7 | `01_Processing/05_Scores/c_data/05_scores.xlsx`, sheet `catalog` (nrow 1110) |
| 06: WGCNA after removing sex and run depth | removeBatchEffect(batch = sex, covariates = log2 run depth) | `06_Modules/README.md` Adjustment; `06_modules.xlsx` read_me (eigengenes "from the sex- and run-depth-adjusted matrix") |
| 06: 827 proteins, 8 modules plus grey | 827 rows; 9 module rows incl. grey (402) | `01_Processing/06_Modules/c_data/06_modules.xlsx`, sheets `modules` (827 rows), `sizes` (turquoise, blue, brown, yellow, green, red, black, pink + grey) |
| "read by 02_Statistics" | proteins.rds, design.rds, scores.rds, modules.rds read by 02_Statistics/01_Fit | step READMEs, "Writes" rows |
| Footer | | `01_Processing/README.md` last paragraph (04_Gene_Sets has no b_reports PDF; the footer is generic) |

Note: Reactome is 140 tested in `04_gene_sets.xlsx` catalog but 139 scored in `05_scores.xlsx` catalog; the image prints no per-collection count.
