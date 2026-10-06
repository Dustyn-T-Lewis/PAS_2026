# 08_Enrichment

Runs every pathway test once, so the arcs only slice the results.

- **Method:** `enrichVolcano::run_enrichment(tests = "fgsea")` on the moderated t of each arc contrast (female PAS, male PAS, treatment × sex), one contrast at a time, seed 20260910. GO:BP and Reactome (disease branch excluded) in the full proteome, MitoCarta in the mito universe, the tested lists from `04_Gene_Sets`; BH within contrast and collection. `enrichVolcano::dedup_terms()` flags redundant terms (EnrichmentMap similarity ≥ 0.375, stronger term kept) for display only. fgsea also runs on GO Slim compartments, on the literature and ROS sets as one list, and on each protein's phenotype-association t per sex. STRING v12 connectivity of the Π and FDR lists.
- **Inputs:** `04_Gene_Sets/c_data/gene_sets.rds`, `07_Fit/c_data/fit.rds`, `07_Fit/c_data/study_full/`, `07_Fit/c_data/study_mito/`
- **Outputs:** `c_data/08_enrichment.xlsx`, `c_data/enrichment.rds` (read by the arcs' 01 to 03), `b_reports/08_enrichment.pdf`
- **Result (p ≤ 0.05 / FDR ≤ 0.05):** female GO:BP 102 / 7, Reactome 17 / 0, MitoCarta 5 / 0; male GO:BP 93 / 4, Reactome 16 / 1, MitoCarta 7 / 3; treatment × sex GO:BP 92 / 0, Reactome 19 / 11, MitoCarta 9 / 6. Tested: 727 GO:BP, 119 Reactome, 43 MitoCarta sets.
- **Note:** each fgsea row carries the same set's singscore contrast from `07_Fit`, tested across mice. fgsea permutes genes, so its FDR runs optimistic against the singscore test. The association fgsea is co-movement across one sex's PAS and vehicle mice, not a PAS effect.
