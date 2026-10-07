# 02_Pathways

Which pathways PAS shifted in males.

| | |
|---|---|
| Reads | `02_Statistics/02_Enrichment/c_data/enrichment.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` |
| Writes | `c_data/02_pathways.xlsx`, `c_data/pathways.rds`, `b_reports/02_pathways.pdf` |
| Run | `Rscript 04_Male/02_Pathways/a_script/02_pathways.R` |

- **Method:** `report_pathways()` keeps `PAS_vs_VEH_male` from `02_Statistics/02_Enrichment`: fgsea on GO:BP and Reactome (full proteome) and MitoCarta (mitochondrial universe), BH within collection, each row beside the same set's singscore contrast. GO Slim BP families organise GO:BP terms for display only.
- **Report:** a volcano ring per collection (12 strongest terms the EnrichmentMap flag keeps), signal by GO Slim family, the 20 strongest full-proteome terms clustered by shared genes, compartments, then the running-score plot of every set at p ≤ 0.05.
- **Result (p ≤ 0.05 / FDR ≤ 0.05):** GO:BP 93 / 4 (lipid homeostasis, plasma membrane organisation, actin filament regulation and supramolecular fibre organisation down); Reactome 16 / 1 (striated muscle contraction down); MitoCarta 7 / 3 (central dogma, translation and mitochondrial ribosome up). The three MitoCarta sets have singscore p ≥ 0.89, so the mouse-level test does not support them.
- **Note:** The fgsea test permutes genes and its FDR runs optimistic; the singscore p beside each row is tested across mice. `fgsea` holds every tested set with its redundancy flag and leading edge; `literature` holds the literature and ROS sets.
