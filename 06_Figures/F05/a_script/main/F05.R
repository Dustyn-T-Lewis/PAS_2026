#!/usr/bin/env Rscript
# Figure 5, what biology carries the sex difference: (a) pathways whose phenotype link differs by
# sex, (b) WGCNA modules, (c) oxidative phosphorylation against soleus pyruvate conductance,
# (d) proteins behind the interaction pathways. Runs the panel scripts in panels/ and assembles
# them.

source(here::here("R", "panels.R"))

panels <- run_panels("F05", "main", c(
  "A_phenotype_links", "B_modules", "C_oxphos", "D_leading_edge"
))
plots <- map(panels, "plot")
tag <- \(p, letter) p + labs(tag = letter)
left <- tag(plots$A_phenotype_links, "a") / tag(plots$C_oxphos, "c") +
  plot_layout(heights = c(1.4, 1))
right <- tag(plots$B_modules, "b") / tag(plots$D_leading_edge, "d") +
  plot_layout(heights = c(1, 1.3))
figure <- (left | right) + plot_layout(widths = c(1.15, 1)) &
  theme(plot.tag = element_text(size = 9, face = "bold"))
save_figure(figure, figure_dirs("F05", "main")$reports, "F05",
  height_mm = 200, width_mm = LANDSCAPE_MM
)
write_panel_data(panels, "F05", "F05")

sessionInfo()
