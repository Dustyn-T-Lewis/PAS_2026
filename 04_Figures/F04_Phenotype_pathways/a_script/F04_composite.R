# Figure 4, phenotype-association pathways: within each sex, the pathway sets that follow each
# phenotype across that sex's mice under the within-treatment model, one per group of overlapping
# sets at fgsea FDR <= 0.05, drawn in the three contrasts and at every trait in both sexes (F03's
# columns). The panel is its own script in panels/; this file runs it and writes one workbook.
# Both models' full association table is in 02_Enrichment's workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- "a_trait_pathways"
panels <- run_panels("F04_Phenotype_pathways", PANELS)

out <- stage_paths("04_Figures", "F04_Phenotype_pathways")
save_figure(
  panels$a_trait_pathways$plot, out$reports, "F04_composite",
  width_mm = 180, height_mm = panels$a_trait_pathways$height_mm
)

write_figure_workbook(
  panels,
  paste(
    "Panel a: every drawn cell: trait block, set, collection, universe, column block, contrast or",
    "sex and trait (within-treatment model), NES, fgsea p, BH FDR, whether that column shows the",
    "set, whether the cell is outlined (the block's trait, shown), the sex with the smaller FDR at",
    "the block's trait (the row's shade), camera p and FDR."
  ),
  file.path(out$data, "F04.xlsx")
)
