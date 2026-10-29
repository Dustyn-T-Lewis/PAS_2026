# Figure 4, phenotype-association pathways: within each sex, the pathway sets that follow each
# phenotype across that sex's mice under the within-treatment model, one per group of overlapping
# sets at fgsea FDR <= 0.05, drawn in the three contrasts and at every trait in both sexes (F03's
# columns), with each trait's distinct significant sets per sex counted beside them. The panel is
# its own script beside this one; this file runs it and writes one workbook, the bar counts its
# own sheet. Both models' full association table is in 02_Enrichment's workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- "a_trait_pathways"
panels <- run_panels("F04_Phenotype_pathways", PANELS)

save_figure(
  panels$a_trait_pathways$plot, "F04_Phenotype_pathways", "F04",
  width_mm = panels$a_trait_pathways$width_mm, height_mm = panels$a_trait_pathways$height_mm
)
panels$counts <- list(data = panels$a_trait_pathways$counts)
write_figure_workbook(panels, "F04_Phenotype_pathways", "F04")
