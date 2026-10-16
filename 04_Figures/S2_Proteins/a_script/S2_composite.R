# Supplementary Figure 2, proteins per contrast: each comparison's Π <= 0.05 proteins across all
# mice, and every protein's fold change and p. Each panel is its own script in panels/; this file
# runs them, stitches them and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("a_pi_female", "b_pi_male", "c_pi_interaction", "d_volcanoes")
panels <- run_panels("S2_Proteins", PANELS)

top <- wrap_plots(map(panels[1:3], "plot"), widths = c(50, 50, 74))
figure <- wrap_elements(full = top) / wrap_elements(full = panels$d_volcanoes$plot) +
  plot_layout(heights = c(130, 62))
out <- stage_paths("04_Figures", "S2_Proteins")
save_figure(figure, out$reports, "S2_composite", width_mm = 180, height_mm = 192)

write_figure_workbook(
  panels,
  c(
    "Panel a: female PAS vs vehicle proteins at Π ≤ 0.05: log2 FC, p, FDR, Π.",
    "Panel b: male PAS vs vehicle proteins at Π ≤ 0.05: log2 FC, p, FDR, Π.",
    "Panel c: treatment × sex proteins at Π ≤ 0.05: log2 FC, p, FDR, Π.",
    "Panel d: every full-proteome protein per contrast: log2 FC, p, FDR, Π, colour class."
  ),
  file.path(out$data, "S2.xlsx")
)
