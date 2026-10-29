# Figure 5, complexes: CORUM complexes and MitoCarta pathways that PAS shifts, against the
# contrasts and the phenotypes, in F03's grid. The panel is its own script beside this one; this
# file runs it and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- "a_set_matrix"
panels <- run_panels("F05_Complexes", PANELS)

save_figure(
  panels$a_set_matrix$plot, "F05_Complexes", "F05",
  width_mm = 180, height_mm = panels$a_set_matrix$height_mm
)
write_figure_workbook(panels, "F05_Complexes", "F05")
