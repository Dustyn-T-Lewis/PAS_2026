# Supplementary Figure 3, proteins per contrast: one column per comparison contrast, every
# protein's fold change and p above that contrast's Π <= 0.05 proteins across all mice and the
# sets each co-expression cluster over-represents. Each panel is its own script beside this one;
# this file runs them, stitches them over one key and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))
source(here("04_Figures", "F02_Proteins", "a_script", "supp", "s3_heatmap.R"))

PANELS <- c("S3a_pi_female", "S3b_pi_male", "S3c_pi_interaction")
panels <- run_panels("F02_Proteins", PANELS)

key <- pi_heatmap_key(unique(list_rbind(map(panels, c("data", "ora")))$collection))
figure <- wrap_plots(map(panels, "plot"), nrow = 1) / wrap_elements(full = key) +
  plot_layout(heights = unit(c(1, 12), c("null", "mm")))
save_figure(figure, "F02_Proteins", "S3", width_mm = 240, height_mm = 186)
write_figure_workbook(panels, "F02_Proteins", "S3")
