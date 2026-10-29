# Supplementary Figure 3, proteins per contrast: each comparison's Π <= 0.05 proteins across all
# mice beside the sets each co-expression cluster over-represents, and every protein's fold change
# and p. Each panel is its own script beside this one; this file runs them, stitches them under
# one key and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("S3a_pi_female", "S3b_pi_male", "S3c_pi_interaction", "S3d_volcanoes")
panels <- run_panels("F02_Proteins", PANELS)

heatmaps <- panels[1:3]
key <- pi_heatmap_key(unique(list_rbind(map(heatmaps, c("data", "ora")))$collection))
figure <- wrap_elements(full = wrap_plots(map(heatmaps, "plot"), nrow = 1)) /
  wrap_elements(full = key) / wrap_elements(full = panels$S3d_volcanoes$plot) +
  plot_layout(heights = c(178, 9, 55))
save_figure(figure, "F02_Proteins", "S3", width_mm = 180, height_mm = 242)
write_figure_workbook(panels, "F02_Proteins", "S3")
