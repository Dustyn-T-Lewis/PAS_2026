# Supplementary Figure 1, processing: what filtering kept, how deep each run went, how the
# normalisation was chosen, what the surrogate variable tracks, the pellet's mitochondrial share,
# and what else went into each pellet. Each panel is its own script beside this one; this file
# runs them, stitches them and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c(
  "S1a_filtering", "S1b_run_depth", "S1c_normalisation", "S1d_surrogate", "S1e_mito_share",
  "S1f_composition"
)
panels <- run_panels("F01_Global", PANELS)

# Panels b, d, e and f share one group key: only f's, at the foot of the figure, is kept. Each row
# lines up its own panels; the rows are stacked whole.
plots <- map(panels, "plot")
keyless <- c("S1b_run_depth", "S1d_surrogate", "S1e_mito_share")
plots[keyless] <- map(plots[keyless], \(p) p + theme(legend.position = "none"))
top <- wrap_plots(plots[1:3], widths = c(60, 45, 50))
bottom <- wrap_plots(plots[4:5], widths = c(80, 70))
figure <- wrap_elements(full = top) / wrap_elements(full = bottom) /
  wrap_elements(full = plots$S1f_composition) +
  plot_layout(heights = c(55, 62, 60))
save_figure(figure, "F01_Global", "S1", width_mm = 180, height_mm = 177)
write_figure_workbook(panels, "F01_Global", "S1")
