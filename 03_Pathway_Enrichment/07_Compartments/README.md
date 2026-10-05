# 07_Compartments

Which cell compartments each contrast shifted, using the GO Slim cellular-component sets.

- **Method:** `enrichVolcano::run_enrichment(tests = c("fgsea", "fry"), min_size = 15, max_size = Inf)`
  on the full proteome, BH within the collection. No upper size limit, so cytosol and nucleus are
  tested. fgsea is competitive; fry asks whether a compartment's proteins move as a whole.
- **Inputs:** `00_Gene_Sets/c_data/gene_sets.rds`,
  `02_Differential_Expression/02_Contrasts/c_data/study_full/da_results.csv`
- **Outputs:** `c_data/07_compartments.xlsx` (sheets: fgsea, fry), `c_data/compartments.rds`
  (fgsea), `c_data/compartments_fry.rds`; no PDFs
- **Result:** 24 of 25 compartments tested. fgsea at FDR ≤ 0.05: female PAS 10 (4 up, 6 down; ER
  up, NES 1.74, FDR 3e-6; nucleoplasm down, NES -1.76); male PAS 2 (cytoskeleton down, NES -1.78;
  ER up); treatment_by_sex 5 (extracellular matrix and region up); sex in VEH 15 and sex in PAS 14,
  with ribosome (NES 3.34 and 3.12) and mitochondrion (NES 1.81 and 1.82) higher in males. fry
  finds 14 per sex contrast and none for treatment.
- **Note:** fgsea reads as a shift in share against the rest of the proteome, not an absolute
  amount. The pellet is mitochondria-enriched and its MitoCarta share differs by sex.
