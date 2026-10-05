#!/usr/bin/env Rscript
# Supplementary Figure 3, compartment shifts with PAS by fry. Runs the panel script in panels/.

source(here::here("R", "panels.R"))

panels <- run_panels("F01", "supp", "S3_A_compartments")
save_figure(panels$S3_A_compartments$plot, figure_dirs("F01", "supp")$reports, "S3",
  height_mm = 100, width_mm = LANDSCAPE_MM
)
write_panel_data(panels, "F01", "S3")

sessionInfo()
