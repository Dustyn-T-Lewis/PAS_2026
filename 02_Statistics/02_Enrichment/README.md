# 02_Enrichment

Runs every pathway test once, so the comparisons only slice the results.

| | |
|---|---|
| Reads | `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `02_Statistics/01_Fit/c_data/fit.rds`, `02_Statistics/01_Fit/c_data/study_full/`, `02_Statistics/01_Fit/c_data/study_mito/` |
| Writes | `c_data/02_enrichment.xlsx`, `c_data/enrichment.rds` (read by `03_Comparisons` and `04_Figures/`), `b_reports/02_enrichment.pdf` |
| Run | `Rscript 02_Statistics/02_Enrichment/a_script/02_enrichment.R` |

- **Contrasts:** `enrichVolcano::run_enrichment(tests = "fgsea")` on the moderated t of each comparison's contrast, one contrast at a time, seed 20260910. GO:BP and Reactome in the full proteome, MitoCarta in the mitochondrial universe; BH within contrast and collection. Each row carries the same set's singscore contrast from `02_Statistics/01_Fit` (log2 FC, p, FDR, r with group), which is tested across mice.
- **Redundancy:** `enrichVolcano::dedup_terms()` flags every term `kept` or `redundant` (EnrichmentMap similarity ≥ 0.375; the stronger term is kept). Figures use the flag to pick display terms; nothing is dropped.
- **Other tests:** GO Slim compartments; the literature and ROS sets as one list, BH across it; fgsea on each protein's phenotype-association t per sex and association model (co-movement, not a PAS effect), beside the singscore association; cross-sex signatures (each sex's p ≤ 0.05 up and down proteins tested in the other sex's ranking); STRING v12 connectivity of the Π and FDR protein lists.
- **Result (p ≤ 0.05 / FDR ≤ 0.05):** female GO:BP 95 / 10, Reactome 24 / 0, MitoCarta 11 / 2; male GO:BP 161 / 35, Reactome 16 / 0, MitoCarta 6 / 0; treatment × sex GO:BP 93 / 20, Reactome 19 / 6, MitoCarta 12 / 11. Three of the four cross-sex signatures carry into the other sex's ranking in the same direction: female up (82 proteins, NES 1.61, p 0.005) and down (95, NES −1.49, p 0.006) in males, and male up in females (74, NES 1.89, p 0.0002); male down does not (116, p 0.52).
- **Basis:** fgsea (Korotkevich et al.) through enrichVolcano for the main tests; the compartment and literature lists call `fgsea::fgsea(nPermSimple = 10000)` directly, as fgsea advises when a set's gene-level statistics are unbalanced, because enrichVolcano keeps the default and such a set would lose its p-value; EnrichmentMap similarity (Merico et al. 2010); Wu & Smyth 2012 on why gene-permutation tests run optimistic when genes are correlated. NES rebuilt from the GSEA definition with 20,000 random sets matches within 0.035.
- **Note:** The fgsea test permutes genes, so its FDR runs optimistic against the singscore test; several FDR-level MitoCarta terms have singscore p > 0.2. The STRING lists are descriptive: subunits of one complex vary together.
