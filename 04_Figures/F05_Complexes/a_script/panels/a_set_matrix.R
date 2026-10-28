# F05 (complexes version) a: the CORUM complexes and MitoCarta pathways that a contrast shows
# (01_Fit's shown flag: singscore p <= 0.05, one set per group of overlapping sets), against the
# three contrasts and against each trait within each sex (pooled model, adjusted for the surrogate
# variable), in F03's grid. Dot fill: r-equivalent of the limma t, so contrasts and traits share
# one scale; area: -log10 p; ring: BH FDR <= 0.05 over every scored set in the universe. Rows
# sit in two blocks by female r-equivalent; labels give the measured members.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(fit = output_path("02_Statistics", "01_Fit", "fit.rds"))
fit <- readRDS(inputs[["fit"]])
BLOCKS <- c(complex = "CORUM", MitoCarta = "MitoCarta")
TONES <- c(SEX_TRAIT_TONES, set_names(names(BLOCKS), BLOCKS))
SHADES <- scales::alpha(c(complex = "grey94", MitoCarta = "grey90"), 0.5)

in_scope <- \(x) filter(x, level == "complex" | collection == "MitoCarta")
contrasts <- in_scope(fit$score_contrasts)
shown <- contrasts |>
  filter(any(shown & contrast %in% COMPARISON_CONTRASTS), .by = feature) |>
  filter(contrast == "PAS_vs_VEH_female") |>
  mutate(key = if_else(level == "complex", "complex", "MitoCarta")) |>
  arrange(match(key, names(BLOCKS)), r_equiv)
rows <- transmute(shown,
  block = factor(BLOCKS[key], BLOCKS), term = feature,
  label = sprintf("%s (%d)", set_label(str_remove(feature, "_\\d+$"), 60), members),
  shade = unname(SHADES[key])
)

cells <- score_cells(contrasts, in_scope(fit$score_associations)) |>
  inner_join(select(rows, feature = term, block, label), by = "feature") |>
  arrange(match(feature, rows$term), section, column) |>
  mutate(column_block = as.character(section), id = paste(section, column))
height_mm <- 3 * nrow(rows) + 32
grid <- pathway_grid(
  transmute(cells,
    term = feature, column = id, value = r_equiv, p = P.Value, selected = adj.P.Val <= ALPHA
  ),
  rows, distinct(cells, label = as.character(column), column = id, column_block),
  title = "Protein complexes and mitochondrial pathways",
  subtitle = sprintf(
    paste(
      "Rows chosen by contrast: singscore p ≤ %s, one set per group of overlapping sets;",
      "measured members in brackets. Traits: pooled association"
    ),
    ALPHA
  ),
  tones = TONES, x_angle = 90, selected_label = paste("FDR ≤", ALPHA),
  fill_name = "r-equivalent", fill_limit = 1, fill_breaks = c(-1, 0, 1), size_name = "p",
  outline_shape = 21
)

save_panel(grid, "F05_Complexes", "a_set_matrix", 180, height_mm)
list(
  plot = grid,
  data = cells |>
    left_join(select(shown, feature, members), by = "feature") |>
    transmute(
      block = as.character(block), feature, label, members, column_block,
      column = as.character(column), r_equiv, P.Value, adj.P.Val
    ),
  inputs = inputs, height_mm = height_mm
)
