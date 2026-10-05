# 01_Scores

One singscore per mouse per gene set, so set-level activity can be classified (03) and related to
phenotypes (04). A score is the mean within-mouse rank of a set's members and tests nothing.

- **Method:** `singscore::rankGenes()` and `multiScore()` on the normalised matrix, one protein per
  gene (the most abundant) with the symbols `enrichVolcano::read_study()` uses. Full-proteome
  collections on the full matrix (sets with ≥ 15 measured members), MitoCarta on the mitochondrial
  matrix (≥ 10).
- **Inputs:** `02_Differential_Expression/02_Contrasts/c_data/contrasts.rds`,
  `02_Differential_Expression/02_Contrasts/c_data/study_full/`, `00_Gene_Sets/c_data/gene_sets.rds`
- **Outputs:** `c_data/01_scores.xlsx` (sheets: catalog, scores_full, scores_mito, plot_index),
  `c_data/set_scores.rds` (read by 03, 04), `b_reports/01_scores_overview.pdf` (one z-scored
  heatmap per collection, display only)
- **Result:** 949 sets scored: 727 GO:BP, 138 Reactome, 34 Hallmark, 7 KEGG, 43 MitoCarta.
