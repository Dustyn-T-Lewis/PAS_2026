# 05_Scores

One singscore per mouse per tested set, complex and cluster: a set's mean within-mouse rank, which uses no group label. 07 tests these scores.

| | |
|---|---|
| Reads | `02_Quantification/c_data/proteins.rds`, `04_Gene_Sets/c_data/gene_sets.rds` |
| Writes | `c_data/05_scores.xlsx`, `c_data/scores.rds` (read by 07), `b_reports/05_scores.pdf` (one z-scored heatmap per collection, display only) |
| Run | `Rscript 00_Trunk/05_Scores/a_script/05_scores.R` |

- **Method:** `singscore::rankGenes()` and `multiScore()`, one protein per gene (the most abundant). Full matrix: GO:BP, Reactome, CORUM, STRING; mitochondrial matrix: MitoCarta.
- **Result:** 995 scores per mouse: 727 GO:BP, 119 Reactome, 93 complexes, 13 clusters, 43 MitoCarta.
