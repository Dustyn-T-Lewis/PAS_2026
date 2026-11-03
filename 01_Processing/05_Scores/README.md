# 05_Scores

One singscore per mouse per tested set, complex and cluster: the mean within-mouse rank of the set's proteins. `02_Statistics/01_Fit` tests these scores.

| | |
|---|---|
| Reads | `01_Processing/02_Quantification/c_data/proteins.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` |
| Writes | `c_data/05_scores.xlsx`, `c_data/scores.rds`, `b_reports/05_scores.pdf` |
| Run | `Rscript 01_Processing/05_Scores/a_script/05_scores.R` |

## Method

- `singscore::rankGenes()` and `multiScore()` (Foroutan et al. 2018), one protein per gene (the most abundant). The score uses no group label.
- Full matrix: GO:BP, Reactome, Hallmark, KEGG, CORUM, STRING. Mitochondrial matrix: MitoCarta.

## Results

- 1,110 scores per mouse: 755 GO:BP, 139 Reactome, 34 Hallmark, 33 KEGG, 98 complexes, 7 clusters, 44 MitoCarta.
