# S8 a: how many proteins each sensitivity variant finds per contrast, against the baseline. Each
# cell holds the count; its colour is the log2 ratio to the baseline (one added to both), so a
# variant that loses results is magenta and one that gains them green.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(sensitivity = output_path("02_Statistics", "07_Sensitivity", "sensitivity.rds"))
counts <- readRDS(inputs[["sensitivity"]])$counts
METRICS <- c(
  p_le_0.05 = "p ≤ 0.05", fdr_le_0.10 = "FDR ≤ 0.10", fdr_le_0.05 = "FDR ≤ 0.05",
  top50_shared = "Top 50 kept"
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
    subtitle = paste(
      "Full proteome, 1,962 proteins (98 expected at p ≤ 0.05 by chance);",
      "top 50 kept: baseline's 50 smallest p still in the variant's 50"
    )
  ) +
  theme_figure() +
  theme(panel.grid.major = element_blank())

save_panel(plot, "S8_Sensitivity", "a_counts", 180, 62)
list(plot = plot, data = counts, inputs = inputs)
