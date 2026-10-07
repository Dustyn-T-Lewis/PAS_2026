# 08_Enrichment

Runs every pathway test once, so the arcs only slice the results.

| | |
|---|---|
| Reads | `04_Gene_Sets/c_data/gene_sets.rds`, `07_Fit/c_data/fit.rds`, `07_Fit/c_data/study_full/`, `07_Fit/c_data/study_mito/` |
| Writes | `c_data/08_enrichment.xlsx`, `c_data/enrichment.rds` (read by the arcs and `Figures/`), `b_reports/08_enrichment.pdf` |
| Run | `Rscript 00_Trunk/08_Enrichment/a_script/08_enrichment.R` |

- **Contrasts:** `enrichVolcano::run_enrichment(tests = "fgsea")` on the moderated t of each arc contrast, one contrast at a time, seed 20260910. GO:BP and Reactome in the full proteome, MitoCarta in the mitochondrial universe; BH within contrast and collection. Each row carries the same set's singscore contrast from 07 (log2 FC, p, FDR, r with group), which is tested across mice.
- **Redundancy:** `enrichVolcano::dedup_terms()` flags every term `kept` or `redundant` (EnrichmentMap similarity ≥ 0.375; the stronger term is kept). Figures use the flag to pick display terms; nothing is dropped.
- **Also:** GO Slim compartments; the literature and ROS sets as one list, BH across it; fgsea on each protein's phenotype-association t per sex (co-movement, not a PAS effect), beside the singscore association; cross-sex signatures (each sex's p ≤ 0.05 up and down proteins tested in the other sex's ranking); STRING v12 connectivity of the Π and FDR protein lists.
- **Result (p ≤ 0.05 / FDR ≤ 0.05):** female GO:BP 102 / 7, Reactome 17 / 0, MitoCarta 5 / 0; male GO:BP 93 / 4, Reactome 16 / 1, MitoCarta 7 / 3; treatment × sex GO:BP 92 / 0, Reactome 19 / 11, MitoCarta 9 / 6. No sex's signature carries into the other's ranking (fgsea p 0.11 to 1).
- **Note:** fgsea permutes genes, so its FDR runs optimistic against the singscore test; several FDR-level MitoCarta terms have singscore p > 0.2. The STRING lists are descriptive: subunits of one complex vary together.
