# Figure 3, contrast pathways: the sets each contrast shows, one per group of overlapping sets at
# fgsea FDR <= 0.05, full proteome under each contrast and MitoCarta last, drawn in all three
# contrasts and at every trait in both sexes (F04's columns), with each contrast's distinct
# significant sets counted beside them. The panel is its own script beside this one; this file
# runs it and writes one workbook, the bar counts its own sheet.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- "a_pathways"
panels <- run_panels("F03_Pathways", PANELS)

save_figure(
  panels$a_pathways$plot, "F03_Pathways", "F03",
  width_mm = panels$a_pathways$width_mm, height_mm = panels$a_pathways$height_mm
)
panels$counts <- list(data = panels$a_pathways$counts)
write_figure_workbook(panels, "F03_Pathways", "F03")
