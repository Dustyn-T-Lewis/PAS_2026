# F01 e: proteins changed per contrast, down to the left and up to the right, the full proteome
# above the mitochondrial universe. Light bar: p or Π <= 0.05; dark bar: Π <= 0.05; dashed: the
# count expected by chance in each direction.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(fit = output_path("02_Statistics", "01_Fit", "fit.rds"))
fit <- readRDS(inputs[["fit"]])
SHORT <- c(PAS_vs_VEH_female = "Female", PAS_vs_VEH_male = "Male", treatment_by_sex = "F − M")

counts <- fit$protein_contrasts |>
  filter(contrast %in% COMPARISON_CONTRASTS) |>
  mutate(direction = if_else(logFC > 0, "Up", "Down")) |>
  summarise(
    tests = n(), nominal = sum(P.Value <= ALPHA | pi_score <= ALPHA), pi = sum(pi_score <= ALPHA),
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
    label = sprintf("%s\nFDR %d", SHORT[contrast], fdr),
    universe = factor(universe_label[universe], universe_label)
  )
chance <- distinct(counts, universe, tests) |>
  mutate(chance = ALPHA / 2 * tests, universe = factor(universe_label[universe], universe_label))
bars <- counts |>
  pivot_longer(c(nominal, pi), names_to = "tier", values_to = "n") |>
  mutate(
    signed = if_else(direction == "Up", n, -n),
    fill = paste0(direction, if_else(tier == "pi", ", Π", ", p or Π")),
    universe = factor(universe_label[universe], universe_label),
    row = factor(contrast, rev(COMPARISON_CONTRASTS))
  ) |>
  arrange(tier == "pi")
fills <- c(
  "Down, p or Π" = DIRECTION_LIGHT[["Down"]], "Up, p or Π" = DIRECTION_LIGHT[["Up"]],
  "Down, Π" = DIRECTION_COLOURS[["Down"]], "Up, Π" = DIRECTION_COLOURS[["Up"]]
)

plot <- ggplot(bars, aes(signed, row)) +
  geom_col(aes(fill = fill), width = 0.7, position = "identity") +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  geom_vline(aes(xintercept = chance), chance, linetype = "dashed", linewidth = 0.25) +
  geom_vline(aes(xintercept = -chance), chance, linetype = "dashed", linewidth = 0.25) +
  contrast_labels(boxes, hjust = 1.06) +
  facet_wrap(~universe, ncol = 1, scales = "free_x") +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = fills, breaks = names(fills), name = NULL) +
  scale_x_continuous(labels = abs, expand = expansion(mult = 0.06)) +
  # Blank tick labels as wide as the boxes reserve their space in the margin.
  scale_y_discrete(labels = \(x) sprintf("%s\nFDR 0", SHORT[x]), expand = expansion(add = 0.5)) +
  labs(
    x = "proteins (down | up)", y = NULL, title = "Proteins per contrast",
    subtitle = "p or Π ≤ 0.05; dashed: expected by chance"
  ) +
  theme_figure() +
  theme(
    panel.grid.major.y = element_blank(), panel.spacing.y = unit(2, "mm"),
    axis.text.y = element_text(
      colour = "transparent", face = "bold", size = 5.9, lineheight = 0.85, margin = margin(r = 3)
    ),
    strip.text = element_text(face = "bold", hjust = 0, margin = margin(b = 1))
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

save_panel(plot, "F01_Global", "e_protein_bars", 60, 62)
list(plot = plot, data = counts, inputs = inputs)
