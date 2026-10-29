# Supplementary Figure 5, sensitivity: how the protein results move when one processing or model
# choice changes. Each panel is its own script beside this one; this file runs them, stitches
# them and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("S5a_counts", "S5b_named")
panels <- run_panels("F02_Proteins", PANELS)

figure <- panels$S5a_counts$plot / panels$S5b_named$plot + plot_layout(heights = c(62, 100))
save_figure(figure, "F02_Proteins", "S5", width_mm = 180, height_mm = 162)
write_figure_workbook(panels, "F02_Proteins", "S5")
