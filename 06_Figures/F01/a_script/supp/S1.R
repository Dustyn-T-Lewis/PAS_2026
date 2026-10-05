#!/usr/bin/env Rscript
# Supplementary Figure 1, quality control: data depth, contamination, normalisation and the
# MitoCarta share per mouse. Runs the panel scripts in panels/ and assembles them.

source(here::here("R", "panels.R"))

panels <- run_panels("F01", "supp", c(
  "S1_A_depth", "S1_B_contamination", "S1_C_normalisation", "S1_D_mito_share"
))
plots <- map(panels, "plot")
figure <- (plots$S1_A_depth + plots$S1_B_contamination) /
  (plots$S1_C_normalisation + plots$S1_D_mito_share + plot_layout(widths = c(1.5, 1))) +
  plot_annotation(tag_levels = "a")
save_figure(figure, figure_dirs("F01", "supp")$reports, "S1",
  height_mm = 150, width_mm = LANDSCAPE_MM
)
write_panel_data(panels, "F01", "S1")

sessionInfo()
