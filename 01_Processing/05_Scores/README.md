# 05_Scores

One singscore per mouse per tested set, complex and cluster: a set's mean within-mouse rank, which uses no group label. `02_Statistics/01_Fit` tests these scores.

| | |
|---|---|
| Reads | `01_Processing/02_Quantification/c_data/proteins.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` |
| Writes | `c_data/05_scores.xlsx`, `c_data/scores.rds` (read by 02_Statistics/01_Fit), `b_reports/05_scores.pdf` (one z-scored heatmap per collection, display only) |
| Run | `Rscript 01_Processing/05_Scores/a_script/05_scores.R` |

- **Method:** `singscore::rankGenes()` and `multiScore()`, one protein per gene (the most abundant). Full matrix: GO:BP, Reactome, CORUM, STRING; mitochondrial matrix: MitoCarta.
- **Result:** 1,043 scores per mouse: 755 GO:BP, 139 Reactome, 98 complexes, 7 clusters, 44 MitoCarta.
- **Basis:** Foroutan et al. 2018 (*BMC Bioinformatics*), singscore's paper: a single-sample score from within-sample ranks, so it needs no group label and no other sample. Rebuilt by hand from that definition, it matches exactly.
