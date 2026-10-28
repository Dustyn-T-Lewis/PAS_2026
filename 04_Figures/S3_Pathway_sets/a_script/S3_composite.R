# Supplementary Figure 3, side lists: the GO Slim compartments and the literature-informed and ROS
# sets by hypothesis, per contrast, in F03's grid. Each panel is its own script in panels/; this
# file runs them, stitches them with one key and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("a_compartments", "b_literature")
panels <- run_panels("S3_Pathway_sets", PANELS)

# One key for both panels (same NES, FDR and database scales). The grids have fixed row heights,
# so patchwork centres the shorter panel b beside a.
figure <- wrap_plots(map(panels, "plot")) + plot_layout(guides = "collect")
out <- stage_paths("04_Figures", "S3_Pathway_sets")
save_figure(figure, out$reports, "S3_composite",
  width_mm = 180, height_mm = panels$a_compartments$height_mm
)

write_figure_workbook(
  panels,
  c(
    "Panel a: GO Slim compartments per contrast: label, NES, p, BH FDR, size.",
    paste(
      "Panel b: literature-informed and ROS sets per contrast: hypothesis, expected direction,",
      "label, database, NES, p, BH FDR over the 19 sets, size."
    )
  ),
  file.path(out$data, "S3.xlsx")
)
