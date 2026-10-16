# Figure 5 (complexes version): CORUM complexes and MitoCarta pathways that PAS shifts, against the
# contrasts and the phenotypes, in the style of the F03 dot matrices. The panel is its own script
# in panels/; this file runs it, checks it against 01_fit.xlsx and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- "a_set_matrix"
panels <- run_panels("F05_Complexes", PANELS)

fit <- readRDS(output_path("02_Statistics", "01_Fit", "fit.rds"))
expected <- fit$score_contrasts |>
  filter(
    contrast %in% COMPARISON_CONTRASTS, level == "complex" | collection == "MitoCarta"
  ) |>
  filter(any(P.Value <= ALPHA), .by = feature) |>
  distinct(feature)
cells <- panels$a_set_matrix$data
stopifnot(
  setequal(cells$feature, expected$feature),
  nrow(cells) == nrow(expected) * (length(COMPARISON_CONTRASTS) + 2 * length(TRAIT_SHORT)),
  !anyNA(cells$r_equiv), !anyNA(cells$members)
)

figure <- panels$a_set_matrix$plot
out <- stage_paths("04_Figures", "F05_Complexes")
save_figure(figure, out$reports, "F05_composite", width_mm = 180, height_mm = 120)

write_figure_workbook(
  panels,
  paste(
    "Panel a: per set and column, the r-equivalent of the limma t, p and BH FDR (over all sets",
    "scored in the universe), with the number of measured members."
  ),
  file.path(out$data, "F05_complexes.xlsx")
)

sessionInfo()
