# Supplementary Figure 4, candidate routes: for each sex, every candidate passing community_routes
# or more routes and which of the six routes it passed. Each panel is its own script beside this
# one; this file runs them, stitches them and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("S4a_female_routes", "S4b_male_routes")
panels <- run_panels("F02_Proteins", PANELS)

figure <- wrap_plots(map(panels, "plot"), ncol = 2)
save_figure(figure, "F02_Proteins", "S4", width_mm = 180, height_mm = 165)
write_figure_workbook(panels, "F02_Proteins", "S4")
