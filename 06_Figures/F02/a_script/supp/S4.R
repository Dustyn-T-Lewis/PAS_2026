#!/usr/bin/env Rscript
# Supplementary Figure 4, behind Figure 2: the PAS-responsive proteins and the literature-informed
# pathways. Runs the panel scripts in panels/ and assembles them.

source(here::here("R", "panels.R"))

panels <- run_panels("F02", "supp", c("S4_A_heatmap", "S4_B_literature"))
plots <- map(panels, "plot")
figure <- wrap_elements(full = plots$S4_A_heatmap) + plots$S4_B_literature +
  plot_layout(widths = c(1.1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))
save_figure(figure, figure_dirs("F02", "supp")$reports, "S4",
  height_mm = 170, width_mm = LANDSCAPE_MM
)
write_panel_data(panels, "F02", "S4")

sessionInfo()
