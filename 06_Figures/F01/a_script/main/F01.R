#!/usr/bin/env Rscript
# Figure 1, the global proteome: (a) PCA with PERMANOVA and the MitoCarta share, (b) baseline sex
# differences, (c) proteins PAS moves, (d) pathways per contrast. Runs the panel scripts in
# panels/ and assembles them.

source(here::here("R", "panels.R"))

panels <- run_panels("F01", "main", c("A_pca", "B_sex_proteins", "C_pas_proteins", "D_pathways"))
plots <- map(panels, "plot")
counts <- plots$B_sex_proteins / plots$C_pas_proteins + plot_layout(heights = c(0.7, 1))
figure <- (plots$A_pca | counts | plots$D_pathways) +
  plot_layout(widths = c(1.1, 1, 1.15)) +
  plot_annotation(tag_levels = "a")
save_figure(figure, figure_dirs("F01", "main")$reports, "F01",
  height_mm = 115,
  width_mm = LANDSCAPE_MM
)
write_panel_data(panels, "F01", "F01")

sessionInfo()
