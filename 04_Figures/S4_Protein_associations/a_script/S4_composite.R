# Supplementary Figure 4, protein associations: how many proteins track each phenotype within each
# sex under each model, and the strongest one per phenotype. Each panel is its own script in
# panels/; this file runs them, checks the counts against 01_fit.xlsx, stitches them and writes one
# workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("a_counts", "b_top_proteins")
panels <- run_panels("S4_Protein_associations", PANELS)

stopifnot(
  nrow(panels$a_counts$data) == length(TRAIT_SHORT) * 2 * 3,
  all(count(panels$b_top_proteins$data, phenotype, sex)$n == length(MODEL_LABELS))
)

figure <- wrap_plots(map(panels, "plot"), ncol = 1, heights = c(60, 70))
out <- stage_paths("04_Figures", "S4_Protein_associations")
save_figure(figure, out$reports, "S4_composite", width_mm = 180, height_mm = 135)

write_figure_workbook(
  panels,
  c(
    "Panel a: per phenotype, sex and model, proteins at p ≤ 0.05, expected by chance, FDR counts.",
    "Panel b: the strongest within-treatment protein per phenotype and sex, under both models."
  ),
  file.path(out$data, "S4.xlsx")
)

sessionInfo()
