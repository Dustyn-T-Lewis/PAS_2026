# Figure 3, contrast pathways: the sets each contrast shows, one per group of overlapping sets at
# fgsea FDR <= 0.05, full proteome under each contrast and MitoCarta last, drawn in all three
# contrasts and at every trait in both sexes (F04's columns), with each contrast's distinct
# significant sets counted beside them. The panel is its own script in panels/; this file runs it
# and writes one workbook, the bar counts its own sheet.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- "a_pathways"
panels <- run_panels("F03_Pathways", PANELS)

out <- stage_paths("04_Figures", "F03_Pathways")
save_figure(
  panels$a_pathways$plot, out$reports, "F03_composite",
  width_mm = panels$a_pathways$width_mm, height_mm = panels$a_pathways$height_mm
)

panels$counts <- list(data = panels$a_pathways$counts, inputs = panels$a_pathways$inputs)
write_figure_workbook(
  panels,
  c(
    paste(
      "Panel a: every drawn cell: row block, set, collection, universe, column block, contrast or",
      "sex and trait (within-treatment model), NES, fgsea p, BH FDR, whether that column shows",
      "the set, whether the cell put the row in its block (outlined), camera p and FDR."
    ),
    paste(
      "Panel a bars: sets at fgsea FDR ≤ 0.05 per block, contrast, direction and collection;",
      "n_distinct (shown, one per group of overlapping sets) is what the bars draw, n_all counts",
      "every set at that FDR. Combinations with no set at FDR ≤ 0.05 are absent."
    )
  ),
  file.path(out$data, "F03.xlsx")
)
