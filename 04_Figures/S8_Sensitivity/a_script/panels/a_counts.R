# S8 a: how many proteins each sensitivity variant finds per contrast, against the baseline. Each
# cell holds the count; its colour is the log2 ratio to the baseline (one added to both), so a
# variant that loses results is magenta and one that gains them green.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(sensitivity = output_path("02_Statistics", "07_Sensitivity", "sensitivity.rds"))
sensitivity <- readRDS(inputs[["sensitivity"]])
counts <- sensitivity$counts
top_n <- sensitivity$top_n
n_proteins <- filter(counts, variant == "baseline", contrast == COMPARISON_CONTRASTS[[1]])$tests
METRICS <- c(
  p_le_0.05 = paste("p ≤", ALPHA), fdr_le_0.10 = "FDR ≤ 0.10", fdr_le_0.05 = paste("FDR ≤", ALPHA),
  top_shared = paste("Top", top_n, "kept")
)

cells <- counts |>
  pivot_longer(all_of(names(METRICS)), names_to = "metric", values_to = "n") |>
  mutate(
    ratio = log2((n + 1) / (n[variant == "baseline"] + 1)),
    .by = c(contrast, metric)
  ) |>
  mutate(
    metric = factor(METRICS[metric], METRICS),
    contrast = factor(contrast, COMPARISON_CONTRASTS),
    variant = factor(variant, rev(levels(variant)))
  )

plot <- ggplot(cells, aes(metric, variant, fill = ratio)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = n), size = 1.9) +
  ggh4x::facet_wrap2(~contrast,
    nrow = 1, labeller = as_labeller(CONTRAST_LABELS),
    strip = contrast_strips(COMPARISON_CONTRASTS)
  ) +
  scale_fill_gradient2(
    low = DIRECTION_COLOURS[["Down"]], high = DIRECTION_COLOURS[["Up"]],
    limits = c(-2.5, 2.5), oob = scales::squish, name = "log2 ratio to baseline",
    guide = guide_colourbar(barwidth = unit(18, "mm"), barheight = unit(2, "mm"))
  ) +
  labs(
    x = NULL, y = NULL, title = "Protein results under each variant",
    subtitle = sprintf(
      "%s proteins, %.0f expected at p ≤ %s; top %d kept: of the baseline's %d smallest p",
      format(n_proteins, big.mark = ","), ALPHA * n_proteins, ALPHA, top_n, top_n
    )
  ) +
  theme_figure() +
  theme(panel.grid.major = element_blank())

save_panel(plot, "S8_Sensitivity", "a_counts", 180, 62)
list(plot = plot, data = counts, inputs = inputs)
