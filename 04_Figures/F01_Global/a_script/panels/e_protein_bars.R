# F01 e: proteins changed per contrast, down to the left and up to the right, the full proteome
# above the mitochondrial universe. Light bar: p <= 0.05; dark bar, drawn over it: Π <= 0.05;
# dashed: the p <= 0.05 count expected by chance in each direction. Each row is named by a box as
# the design panel draws it: the contrast's tint inside a border in its colour, with its FDR count.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(fit = output_path("02_Statistics", "01_Fit", "fit.rds"))
fit <- readRDS(inputs[["fit"]])

counts <- fit$protein_contrasts |>
  filter(contrast %in% COMPARISON_CONTRASTS) |>
  mutate(direction = if_else(logFC > 0, "Up", "Down")) |>
  summarise(
    tests = n(), nominal = sum(P.Value <= ALPHA), pi = sum(pi_score <= ALPHA),
    fdr = sum(adj.P.Val <= ALPHA), .by = c(universe, contrast, direction)
  ) |>
  mutate(tests = sum(tests), .by = c(universe, contrast))
universe_label <- distinct(counts, universe, tests) |>
  mutate(label = sprintf(
    "%s (%s proteins)", c(full = "Full proteome", mito = "Mitochondrial")[universe],
    format(tests, big.mark = ",", trim = TRUE)
  )) |>
  select(universe, label) |>
  deframe()
boxes <- counts |>
  summarise(fdr = sum(fdr), .by = c(universe, contrast)) |>
  transmute(
    x = -Inf, y = factor(contrast, rev(COMPARISON_CONTRASTS)), contrast,
    label = sprintf("%s\nFDR %d", CONTRAST_SHORT[contrast], fdr),
    universe = factor(universe_label[universe], universe_label)
  )
chance <- distinct(counts, universe, tests) |>
  mutate(chance = ALPHA / 2 * tests, universe = factor(universe_label[universe], universe_label))
tiers <- c(nominal = paste("p ≤", ALPHA), pi = paste("Π ≤", ALPHA))
bars <- counts |>
  pivot_longer(c(nominal, pi), names_to = "tier", values_to = "n") |>
  mutate(
    signed = if_else(direction == "Up", n, -n),
    fill = paste0(direction, ", ", tiers[tier]),
    universe = factor(universe_label[universe], universe_label),
    row = factor(contrast, rev(COMPARISON_CONTRASTS))
  ) |>
  arrange(tier == "pi")
fills <- c(
  set_names(DIRECTION_LIGHT, paste0(names(DIRECTION_LIGHT), ", ", tiers[["nominal"]])),
  set_names(DIRECTION_COLOURS, paste0(names(DIRECTION_COLOURS), ", ", tiers[["pi"]]))
)

plot <- ggplot(bars, aes(signed, row)) +
  geom_col(aes(fill = fill), width = 0.7, position = "identity") +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  geom_vline(aes(xintercept = chance), chance, linetype = "dashed", linewidth = 0.25) +
  geom_vline(aes(xintercept = -chance), chance, linetype = "dashed", linewidth = 0.25) +
  geom_label(
    aes(x, y,
      label = label, fill = I(contrast_tone("fill")[contrast]),
      colour = I(contrast_tone("line")[contrast])
    ),
    boxes,
    hjust = 1.06, text.colour = "black", fontface = "bold", size = 1.8, linewidth = 0.4,
    label.r = unit(0.5, "mm"), label.padding = unit(0.45, "mm"), lineheight = 0.85,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  facet_wrap(~universe, ncol = 1, scales = "free_x") +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = fills, breaks = names(fills), name = NULL) +
  scale_x_continuous(labels = abs, expand = expansion(mult = 0.06)) +
  # Blank tick labels as wide as the boxes reserve their space in the margin.
  scale_y_discrete(
    labels = \(x) sprintf("%s\nFDR 0", CONTRAST_SHORT[x]), expand = expansion(add = 0.5)
  ) +
  labs(
    x = "proteins (down / up)", y = NULL, title = "Proteins per contrast",
    subtitle = sprintf("Dashed: p ≤ %s expected by chance", ALPHA)
  ) +
  theme_figure() +
  theme(
    panel.grid.major.y = element_blank(), panel.spacing.y = unit(2, "mm"),
    axis.text.y = element_text(
      colour = "transparent", face = "bold", size = 5.9, lineheight = 0.85, margin = margin(r = 3)
    ),
    strip.text = element_text(hjust = 0, margin = margin(b = 1))
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

save_panel(plot, "F01_Global", "e_protein_bars", 57, 64)
list(plot = plot, data = counts, inputs = inputs)
