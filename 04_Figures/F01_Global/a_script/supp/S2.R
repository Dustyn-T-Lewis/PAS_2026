# Supplementary Figure 2, structure: what the main proteome components track, every WGCNA module
# against the contrasts and the phenotypes, and whether more features separate PAS from vehicle
# perfectly than under relabelling of the mice. Each panel is its own script beside this one; this
# file runs them, stitches them and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("S2a_component_r2", "S2b_module_traits", "S2c_perfect_separators")
panels <- run_panels("F01_Global", PANELS)

# The two half-width panels share the top row; the module heatmap takes the full width below.
top <- wrap_plots(panels$S2a_component_r2$plot, panels$S2c_perfect_separators$plot)
figure <- wrap_elements(full = top) / panels$S2b_module_traits$plot +
  plot_layout(heights = c(60, 92))
save_figure(figure, "F01_Global", "S2", width_mm = 180, height_mm = 152)
write_figure_workbook(panels, "F01_Global", "S2")
