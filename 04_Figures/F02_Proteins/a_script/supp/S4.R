# Supplementary Figure 4, sensitivity: how the protein results move when one processing or model
# choice changes. Each panel is its own script beside this one; this file runs them, stitches
# them and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("S4a_counts", "S4b_named")
panels <- run_panels("F02_Proteins", PANELS)

figure <- panels$S4a_counts$plot / panels$S4b_named$plot + plot_layout(heights = c(62, 100))
save_figure(figure, "F02_Proteins", "S4", width_mm = 180, height_mm = 162)
write_figure_workbook(panels, "F02_Proteins", "S4")
