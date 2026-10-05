#!/usr/bin/env Rscript
# Figure 2, which pathways PAS moves in each sex: (a) full-proteome rings, (b) mitochondrial rings,
# female above male, (c) the term tree with each sex and the interaction. Runs the panel scripts in
# panels/ and assembles them.

source(here::here("R", "panels.R"))

panels <- run_panels("F02", "main", c("A_full_rings", "B_mito_rings", "C_tree"))
plots <- map(panels, "plot")
figure <- wrap_elements(full = plots$A_full_rings) + wrap_elements(full = plots$B_mito_rings) +
  plots$C_tree + plot_layout(widths = c(1, 1, 1.25)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))
save_figure(figure, figure_dirs("F02", "main")$reports, "F02",
  height_mm = 175, width_mm = LANDSCAPE_MM
)
write_panel_data(panels, "F02", "F02")

sessionInfo()
