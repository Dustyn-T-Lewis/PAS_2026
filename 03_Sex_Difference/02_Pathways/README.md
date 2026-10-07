# 02_Pathways

Which pathways respond to PAS differently in females and males.

| | |
|---|---|
| Reads | `00_Trunk/08_Enrichment/c_data/enrichment.rds`, `00_Trunk/04_Gene_Sets/c_data/gene_sets.rds` |
| Writes | `c_data/02_pathways.xlsx`, `c_data/pathways.rds`, `b_reports/02_pathways.pdf` |
| Run | `Rscript 03_Sex_Difference/02_Pathways/a_script/02_pathways.R` |

- **Method:** `arc_pathways()` keeps `treatment_by_sex` from `08_Enrichment`. A positive NES means PAS raised the set more in females than in males.
- **Report:** as in the female and male arcs: rings, GO Slim families, the term tree, compartments, then the running-score plot of every set at p ≤ 0.05.
- **Result (p ≤ 0.05 / FDR ≤ 0.05):** GO:BP 92 / 0; Reactome 19 / 11 (extracellular matrix, collagen formation, ECM proteoglycans, scavenger receptors, chemical stress and striated muscle contraction higher in the female response; cell cycle lower); MitoCarta 9 / 6 (central dogma, mitochondrial ribosome and translation lower in the female response; metabolism, carbohydrate and lipid metabolism higher).
- **Note:** the six MitoCarta sets have singscore p 0.14 to 0.54, so the mouse-level test does not support them; seven of the 11 Reactome sets have singscore p ≤ 0.05 (cell cycle, mitotic cell cycle, collagen formation, ECM proteoglycans, ECM degradation, muscle contraction, striated muscle contraction).
