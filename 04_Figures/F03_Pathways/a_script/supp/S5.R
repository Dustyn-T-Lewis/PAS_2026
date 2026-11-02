# Supplementary Figure 5, side lists: the GO Slim compartments and the literature-informed and ROS
# sets by hypothesis, per contrast, in F03's grid. Each panel is its own script beside this one;
# this file runs them, stitches them with one key and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("S5a_compartments", "S5b_literature")
panels <- run_panels("F03_Pathways", PANELS)

# One key for both panels (same NES, FDR and database scales). The grids have fixed row heights,
# so patchwork centres the shorter panel b beside a.
figure <- wrap_plots(map(panels, "plot")) + plot_layout(guides = "collect")
save_figure(figure, "F03_Pathways", "S5",
  width_mm = 180, height_mm = panels$S5a_compartments$height_mm
)
write_figure_workbook(panels, "F03_Pathways", "S5")
