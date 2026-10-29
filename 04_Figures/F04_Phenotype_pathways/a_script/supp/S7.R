# Supplementary Figure 7, protein associations: how many proteins track each phenotype within each
# sex under each model and in which direction, and the three strongest per phenotype. Each panel
# is its own script beside this one; this file runs them, stitches them and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("S7a_counts", "S7b_top_proteins")
panels <- run_panels("F04_Phenotype_pathways", PANELS)

# Panel a's axis title sits by its own axis rather than level with b's longer labels.
figure <- wrap_plots(
  free(panels$S7a_counts$plot, type = "label", side = "l"), panels$S7b_top_proteins$plot,
  ncol = 1, heights = c(62, 100)
)
save_figure(figure, "F04_Phenotype_pathways", "S7", width_mm = 180, height_mm = 165)
write_figure_workbook(panels, "F04_Phenotype_pathways", "S7")
