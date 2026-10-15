# Helpers the four steps in 03_Comparisons/code/ share: where a comparison writes, which mice it
# looks at, its colours and the report grids. Sourced by every step; it sources helpers.R.

source(here::here("00_Setup", "Shared_functions", "helpers.R"))

GRID_PER_PAGE <- 20
LEVEL_LABELS <- c(
  protein = "Proteins", pathway = "Pathways", module = "Modules", complex = "Complexes",
  cluster = "Clusters"
)

# A comparison writes every step into one folder: b_reports/ and c_data/, files named by step.
comparison_paths <- function(comparison, step) {
  root <- here("03_Comparisons", COMPARISONS[[comparison]]$dir)
  walk(file.path(root, c("b_reports", "c_data")), dir.create, showWarnings = FALSE)
  list(
    root = root, data = file.path(root, "c_data"), reports = file.path(root, "b_reports"),
    workbook = file.path(root, "c_data", paste0(step, ".xlsx"))
  )
}

# The mice a comparison looks at: one sex for the PAS comparisons, all twenty for the sex
# difference.
comparison_mice <- function(comparison, targets) {
  sex_code <- COMPARISONS[[comparison]]$sex
  if (is.na(sex_code)) targets$sample_id else targets$sample_id[targets$sex == sex_code]
}

comparison_colour <- function(comparison) {
  unname(CONTRAST_COLOURS[CONTRAST_LABELS[COMPARISONS[[comparison]]$contrast]])
}

# The strongest tier a protein reaches: FDR, then Π, then nominal p.
tier_of <- function(p, fdr, pi) {
  factor(
    case_when(
      fdr <= ALPHA ~ "FDR ≤ 0.05", pi <= ALPHA ~ "Π ≤ 0.05",
      p <= ALPHA ~ "p ≤ 0.05",
      .default = "p > 0.05"
    ),
    c("p > 0.05", "p ≤ 0.05", "Π ≤ 0.05", "FDR ≤ 0.05")
  )
}
tier_colours <- function(comparison) {
  colour <- comparison_colour(comparison)
  c(
    "p > 0.05" = "grey82", "p ≤ 0.05" = lighten(colour, 0.55),
    "Π ≤ 0.05" = colour, "FDR ≤ 0.05" = "black"
  )
}

# Splits long data into report pages of per_page facets, keeping the order of the strip factor; a
# short last page is padded with blank facets so every page keeps the same grid.
paginate <- function(data, draw, per_page = GRID_PER_PAGE) {
  strips <- levels(data$strip)
  pad <- strrep(" ", seq_len(per_page))
  split(strips, ceiling(seq_along(strips) / per_page)) |>
    map(\(s) {
      facet_levels <- c(s, pad[seq_len(per_page - length(s))])
      draw(filter(data, strip %in% s) |> mutate(strip = factor(strip, facet_levels)))
    })
}

# Facet labels in first-seen order.
strip_factor <- function(label) factor(label, unique(label))
short <- \(x, width = 34) str_trunc(enrichVolcano::clean_label(x, width = 200), width)

# Values by mouse for many features, one facet each, grouped and coloured by group.
mouse_grid <- function(long, title, subtitle, y = "log2 abundance") {
  paginate(long, \(d) {
    ggplot(d, aes(group, value, colour = group, shape = sex)) +
      stat_summary(
        fun = median, geom = "crossbar", width = 0.5, linewidth = 0.25, colour = "grey40"
      ) +
      geom_point(position = position_jitter(width = 0.12, seed = SEED), size = 1.1) +
      facet_wrap(~strip, scales = "free_y", ncol = 5, nrow = GRID_PER_PAGE / 5, drop = FALSE) +
      scale_colour_manual(values = GROUP_COLOURS, guide = "none") +
      scale_shape_manual(values = SEX_SHAPES, guide = "none") +
      labs(x = NULL, y = y, title = title, subtitle = subtitle) +
      theme_figure() +
      theme(strip.text = element_text(size = 5, face = "plain"), axis.text = element_text(size = 5))
  })
}
long_values <- function(m, features, strips, targets) {
  as_tibble(m[features, targets$sample_id, drop = FALSE], rownames = "feature") |>
    mutate(strip = strips) |>
    pivot_longer(-c(feature, strip), names_to = "sample_id") |>
    left_join(select(targets, sample_id, group, sex), by = "sample_id") |>
    mutate(group = factor(group, FIGURE_GROUPS), strip = strip_factor(strip))
}
