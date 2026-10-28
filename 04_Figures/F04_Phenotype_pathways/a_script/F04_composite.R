# Figure 4, phenotype-association pathways: within each sex, the pathway sets that follow each
# phenotype across that sex's mice under the within-treatment model, one per group of overlapping
# sets at fgsea FDR <= 0.05, drawn in the three contrasts and at every trait in both sexes (F03's
# columns), with each trait's distinct significant sets per sex counted beside them. The panel is
# its own script in panels/; this file runs it and writes one workbook, the bar counts its own
# sheet. Both models' full association table is in 02_Enrichment's workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- "a_trait_pathways"
panels <- run_panels("F04_Phenotype_pathways", PANELS)

out <- stage_paths("04_Figures", "F04_Phenotype_pathways")
save_figure(
  panels$a_trait_pathways$plot, out$reports, "F04_composite",
  width_mm = panels$a_trait_pathways$width_mm, height_mm = panels$a_trait_pathways$height_mm
)

panels$counts <- list(
  data = panels$a_trait_pathways$counts, inputs = panels$a_trait_pathways$inputs
)
write_figure_workbook(
  panels,
  c(
    paste(
      "Panel a: every drawn cell: trait block, set, collection, universe, column block, contrast",
      "or sex and trait (within-treatment model), NES, fgsea p, BH FDR, whether that column shows",
      "the set, whether the cell is outlined (the block's trait, shown), the sex that chose the",
      "row (its shade), camera p and FDR."
    ),
    paste(
      "Panel a bars: within-treatment sets at fgsea FDR ≤ 0.05 per trait, sex, direction (NES",
      "sign) and collection; n_distinct (shown, one per group of overlapping sets) is what the",
      "bars draw, n_all counts every set at that FDR. Combinations with no set at FDR ≤ 0.05 are",
      "absent."
    )
  ),
  file.path(out$data, "F04.xlsx")
)
