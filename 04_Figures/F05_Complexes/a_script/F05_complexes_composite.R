# Figure 5 (complexes version): CORUM complexes and MitoCarta pathways that PAS shifts, against the
# contrasts and the phenotypes, in the style of the F03 dot matrices. The panel is its own script
# in panels/; this file runs it and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- "a_set_matrix"
panels <- run_panels("F05_Complexes", PANELS)

out <- stage_paths("04_Figures", "F05_Complexes")
save_figure(
  panels$a_set_matrix$plot, out$reports, "F05_complexes_composite",
  width_mm = 180, height_mm = 120
)

write_figure_workbook(
  panels,
  paste(
    "Panel a: per set and column, the r-equivalent of the limma t, p and BH FDR (over all sets",
    "scored in the universe), with the number of measured members."
  ),
  file.path(out$data, "F05_complexes.xlsx")
)
