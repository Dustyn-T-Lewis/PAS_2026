# Figure 3, contrast pathways: the sets each contrast shows, one per group of overlapping sets at
# fgsea FDR <= 0.05, full proteome under each contrast and MitoCarta last, drawn in all three
# contrasts and at every trait in both sexes (F04's columns). The panel is its own script in
# panels/; this file runs it and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- "a_pathways"
panels <- run_panels("F03_Pathways", PANELS)

out <- stage_paths("04_Figures", "F03_Pathways")
save_figure(
  panels$a_pathways$plot, out$reports, "F03_composite",
  width_mm = 180, height_mm = panels$a_pathways$height_mm
)

write_figure_workbook(
  panels,
  paste(
    "Panel a: every drawn cell: row block, set, collection, universe, column block, contrast or",
    "sex and trait (within-treatment model), NES, fgsea p, BH FDR, whether that column shows the",
    "set, whether the cell put the row in its block (outlined), camera p and FDR."
  ),
  file.path(out$data, "F03.xlsx")
)
