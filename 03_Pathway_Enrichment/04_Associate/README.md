# 04_Associate

Which gene sets track each of the six phenotypes, read as NES, FDR and leading edge like the
treatment results.

- **Method:** Lead: `enrichVolcano::run_enrichment(tests = "fgsea", min_size = 15)` on proteins
  ranked by their phenotype t, then `dedup_terms(method = "enrichmentmap")`. Supplementary: `limma`
  on the singscores, BH within phenotype, model and collection. Three models from
  `02_Differential_Expression/04_Associate`: pooled slope, females only, slope difference (male
  minus female).
- **Inputs:** `02_Differential_Expression/04_Associate/c_data/protein_association.rds`,
  `00_Gene_Sets/c_data/gene_sets.rds`, `01_Scores/c_data/set_scores.rds`,
  `02_Differential_Expression/02_Contrasts/c_data/contrasts.rds`, `00_Input/phenotypes.xlsx`
- **Outputs:** `c_data/04_associate.xlsx` (sheets: summary, phenotype_gsea, singscore,
  plot_index), `b_reports/04_associate_<phenotype>.pdf` (12 smallest-p terms per model)
- **Result:** fgsea gives 1,223 phenotype-model-term hits at FDR ≤ 0.05 across 18 phenotype-model
  runs of 939 sets. Largest: conductance_sol_pyr slope difference 271 (223 up), conductance_sol_pyr
  females only 205, capillary_density pooled 176, conductance_pla_pyr pooled 106. The singscore
  limma check finds 1 (a Hallmark set, conductance_sol_oct, females only).
- **Note:** The rankings come from all 1,842 proteins, so MitoCarta is tested against the full
  proteome with `min_size` 15 (33 sets), not in the mitochondrial universe as in 02. fgsea and the
  singscore check disagree by three orders of magnitude; treat the fgsea counts as a ranking.
