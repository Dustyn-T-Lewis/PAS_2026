# 05_Volcano

Draws each contrast's protein volcano with its strongest fgsea terms around it, one collection per
page. Display only.

- **Method:** `enrichVolcano::plot_volcano_ring(collapse = TRUE, term_threshold = 0.05, n_terms = 12)`:
  up to 12 terms at FDR ≤ 0.05, one per EnrichmentMap cluster. Pages: Hallmark, Reactome, KEGG,
  GO:BP (full proteome), MitoCarta (mitochondrial universe).
- **Inputs:** `02_Contrasts/c_data/enrichment.rds`,
  `02_Differential_Expression/02_Contrasts/c_data/study_full/` and `study_mito/`
- **Outputs:** `c_data/05_volcano.xlsx` (sheet: plot_index),
  `b_reports/05_volcano_<contrast>.pdf` (five, 8 × 8 in)
- **Result:** No new statistics; term counts are those in `02_Contrasts`.
