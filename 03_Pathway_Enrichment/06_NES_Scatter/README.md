# 06_NES_Scatter

Plots pathway effects in one contrast against another: female PAS against male PAS, and sex in VEH
against sex in PAS. Display only.

- **Method:** `enrichVolcano::plot_scatter(comparison = "concordance", collapse = TRUE, term_threshold = 0.05, label_n = 10)`
  on fgsea NES. Pages: Hallmark, Reactome, GO:BP, MitoCarta (no KEGG).
- **Inputs:** `02_Contrasts/c_data/enrichment.rds`
- **Outputs:** `c_data/06_nes_scatter.xlsx` (sheet: plot_index),
  `b_reports/06_nes_scatter_treatment.pdf`, `b_reports/06_nes_scatter_sex.pdf`
- **Result:** No new statistics.
- **Note:** Whether the sexes respond differently is tested by `treatment_by_sex` in
  `02_Contrasts`, not read off these plots.
