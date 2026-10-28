# Supplementary Figure 1, processing: what filtering kept, how deep each run went, how the
# normalisation was chosen, what the surrogate variable tracks, the pellet's mitochondrial share,
# and what else went into each pellet.
# Each panel is its own script in panels/; this file runs them, stitches them and writes one
# workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c(
  "a_filtering", "b_run_depth", "c_normalisation", "d_surrogate", "e_mito_share", "f_composition"
)
panels <- run_panels("S1_Processing", PANELS)

# Panels b, d, e and f share one group key: only f's, at the foot of the figure, is kept. Each row
# lines up its own panels; the rows are stacked whole.
plots <- map(panels, "plot")
keyless <- c("b_run_depth", "d_surrogate", "e_mito_share")
plots[keyless] <- map(plots[keyless], \(p) p + theme(legend.position = "none"))
top <- wrap_plots(plots[1:3], widths = c(60, 45, 50))
bottom <- wrap_plots(plots[4:5], widths = c(80, 70))
figure <- wrap_elements(full = top) / wrap_elements(full = bottom) /
  wrap_elements(full = plots$f_composition) +
  plot_layout(heights = c(55, 62, 60))
out <- stage_paths("04_Figures", "S1_Processing")
save_figure(figure, out$reports, "S1_composite", width_mm = 180, height_mm = 177)

write_figure_workbook(
  panels,
  c(
    "Panel a: precursors and protein groups entering filtering and remaining after each filter.",
    "Panel b: precursors and proteins identified per run.",
    "Panel c: PCV, PMAD and PEV per normalisation method and universe; eligibility; choice.",
    "Panel d: surrogate variable 1 and log2 run depth per mouse.",
    "Panel e: % of signal from MitoCarta proteins per mouse.",
    "Panel f: per mouse and marker group shown, the composition score."
  ),
  file.path(out$data, "S1.xlsx")
)
