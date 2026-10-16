# Supplementary Figure 4, protein associations: how many proteins track each phenotype within each
# sex under each model, and the strongest one per phenotype. Each panel is its own script in
# panels/; this file runs them, stitches them and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("a_counts", "b_top_proteins")
panels <- run_panels("S4_Protein_associations", PANELS)

# Panel a's axis title sits by its own axis rather than level with b's longer labels.
figure <- wrap_plots(
  free(panels$a_counts$plot, type = "label", side = "l"), panels$b_top_proteins$plot,
  ncol = 1, heights = c(60, 70)
)
out <- stage_paths("04_Figures", "S4_Protein_associations")
save_figure(figure, out$reports, "S4_composite", width_mm = 180, height_mm = 135)

write_figure_workbook(
  panels,
  c(
    "Panel a: per phenotype, sex and model, proteins at p ≤ 0.05, expected by chance, FDR counts.",
    "Panel b: the strongest within-treatment protein per phenotype and sex, under both models."
  ),
  file.path(out$data, "S4.xlsx")
)
