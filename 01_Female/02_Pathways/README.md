# 02_Pathways

Which pathways PAS shifted in females.

| | |
|---|---|
| Reads | `00_Trunk/08_Enrichment/c_data/enrichment.rds`, `00_Trunk/04_Gene_Sets/c_data/gene_sets.rds` |
| Writes | `c_data/02_pathways.xlsx`, `c_data/pathways.rds`, `b_reports/02_pathways.pdf` |
| Run | `Rscript 01_Female/02_Pathways/a_script/02_pathways.R` |

- **Method:** `arc_pathways()` keeps `PAS_vs_VEH_female` from `08_Enrichment`: fgsea on GO:BP and Reactome (full proteome) and MitoCarta (mitochondrial universe), BH within collection, each row beside the same set's singscore contrast. GO Slim BP families organise GO:BP terms for display only.
- **Report:** a volcano ring per collection (12 strongest terms the EnrichmentMap flag keeps), signal by GO Slim family, the 20 strongest full-proteome terms clustered by shared genes, compartments, then the running-score plot of every set at p ≤ 0.05.
- **Result (p ≤ 0.05 / FDR ≤ 0.05):** GO:BP 102 sets at p ≤ 0.05 (36 expected), 7 at FDR ≤ 0.05 (muscle system process, myofibril assembly, sarcomere organisation and carbohydrate catabolism up; cell cycle down); Reactome 17 / 0; MitoCarta 5 / 0. None of the 7 FDR-level GO:BP sets is matched by its singscore contrast except cell cycle (p 0.039).
- **Note:** fgsea permutes genes and its FDR runs optimistic; the singscore p beside each row is tested across mice. `fgsea` holds every tested set with its redundancy flag and leading edge; `literature` holds the literature and ROS sets.
