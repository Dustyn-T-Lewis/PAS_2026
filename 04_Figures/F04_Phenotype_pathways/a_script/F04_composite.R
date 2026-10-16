# Figure 4, phenotype-association pathways: within each sex, the pathway sets that follow each
# phenotype across that sex's mice under the within-treatment model, one per group of overlapping
# sets at fgsea FDR <= 0.05. The panel is its own script in panels/; this file runs it and writes
# one workbook. Both models' full association table is in 02_Enrichment's workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- "a_trait_pathways"
panels <- run_panels("F04_Phenotype_pathways", PANELS)

out <- stage_paths("04_Figures", "F04_Phenotype_pathways")
# Wrapped like the multi-panel composites, so the letter keeps the same margin from the edge.
save_figure(
  wrap_plots(panels$a_trait_pathways$plot), out$reports, "F04_composite",
  width_mm = 180, height_mm = panels$a_trait_pathways$height_mm
)

write_figure_workbook(
  panels,
  paste(
    "Panel a: every drawn cell (within-treatment model): sex, trait block, set, collection,",
    "universe, phenotype, NES, fgsea p and BH FDR, whether that phenotype shows the set, camera",
    "p and FDR."
  ),
  file.path(out$data, "F04.xlsx")
)
