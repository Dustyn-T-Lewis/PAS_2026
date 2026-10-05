# 02_Contrasts

Tests which pathways each of the five contrasts moved. fgsea is the lead result; camera
(competitive) and fry (self-contained) check it.

- **Method:** `enrichVolcano::run_enrichment(tests = c("fgsea", "camera", "fry"), inter_gene_cor = 0.01, max_size = 500)`
  on the limpa results, precision weights, one protein per gene (the most abundant), BH within
  collection, seed fixed. `min_size` 15 for the full proteome, 10 for MitoCarta. Symbols from
  `org.Mm.eg.db` (1,804 of 1,842 proteins map; the rest keep DIA-NN's). `dedup_terms(method =
  "enrichmentmap")` adds a redundancy flag (`dedup_status`) for display.
- **Inputs:** `00_Gene_Sets/c_data/gene_sets.rds`,
  `02_Differential_Expression/02_Contrasts/c_data/study_full/` and `study_mito/` (`da_results.csv`)
- **Outputs:** `c_data/02_contrasts.xlsx` (sheets: summary, fgsea, camera, fry, ros_sets,
  plot_index), `c_data/enrichment.rds` (read by 05, 06, 08),
  `b_reports/02_contrasts_<contrast>.pdf` (15 smallest-p fgsea terms per collection)
- **Result (fgsea, FDR ≤ 0.05):** female PAS 15 terms (8 up, 7 down: 7 GO:BP, 5 KEGG, 2 Reactome,
  1 Hallmark); male PAS 9 (3 up, 6 down: 4 GO:BP, 3 MitoCarta, 2 Hallmark); treatment_by_sex 18
  (13 up, 5 down: 10 Reactome, 6 MitoCarta, 2 Hallmark); sex_in_VEH 69; sex_in_PAS 68. fry finds
  5 Reactome terms (all down) in female PAS and none in male PAS or the interaction.
- **Result (ROS sets):** none passes FDR. MitoCarta ROS and glutathione metabolism rises in female
  PAS (NES 1.73, p = 0.013, FDR 0.16); Reactome ROS detoxification has too few measured proteins
  to test.
- **Note:** Tested per collection: 727 GO:BP, 138 Reactome, 34 Hallmark, 43 MitoCarta, and only 7
  of 642 KEGG MEDICUS sets reach 15 measured members.
