# F01 e: proteins changed per contrast, down to the left and up to the right, in the full proteome
# and the mitochondrial universe. Light bar: p or Π <= 0.05; dark bar: Π <= 0.05; dashed: chance.
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
universe_size <- deframe(distinct(counts, universe, tests))
universe_label <- c(
  full = sprintf("Full proteome (%s)", format(universe_size[["full"]], big.mark = ",")),
  mito = sprintf("Mitochondrial (%s)", format(universe_size[["mito"]], big.mark = ","))
)
row_label <- counts |>
  filter(universe == "full") |>
  summarise(fdr = sum(fdr), .by = contrast) |>
  mutate(label = sprintf("%s\nFDR %d", SHORT[contrast], fdr)) |>
  select(contrast, label) |>
  deframe()
row_label <- row_label[COMPARISON_CONTRASTS]
chance <- distinct(counts, universe, tests) |>
  mutate(chance = ALPHA / 2 * tests, universe = universe_label[universe])
bars <- counts |>
  pivot_longer(c(nominal, pi), names_to = "tier", values_to = "n") |>
  mutate(
    signed = if_else(direction == "Up", n, -n),
    fill = paste0(direction, if_else(tier == "pi", ", Π", ", p or Π")),
    universe = factor(universe_label[universe], universe_label),
    row = factor(row_label[contrast], rev(row_label))
  ) |>
  arrange(tier == "pi")
boxes <- tibble(
  x = -Inf, y = factor(row_label, rev(row_label)), label = row_label, contrast = names(row_label),
  universe = factor(universe_label[["full"]], universe_label)
)
fills <- c(
  "Down, p or Π" = DIRECTION_LIGHT[["Down"]], "Down, Π" = DIRECTION_COLOURS[["Down"]],
  "Up, p or Π" = DIRECTION_LIGHT[["Up"]], "Up, Π" = DIRECTION_COLOURS[["Up"]]
)

plot <- ggplot(bars, aes(signed, row)) +
  geom_col(aes(fill = fill), width = 0.7, position = "identity") +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  geom_vline(aes(xintercept = chance), chance, linetype = "dashed", linewidth = 0.25) +
  geom_vline(aes(xintercept = -chance), chance, linetype = "dashed", linewidth = 0.25) +
  contrast_labels(boxes, hjust = 1.04) +
  facet_wrap(~universe, scales = "free_x") +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = fills, breaks = names(fills), name = NULL) +
  scale_x_continuous(labels = abs, expand = expansion(mult = 0.08)) +
  labs(
    x = "proteins (down | up)", y = NULL, title = "Proteins per contrast",
    subtitle = "Light: p or Π ≤ 0.05; dark: Π ≤ 0.05; dashed: chance"
  ) +
  theme_figure() +
  theme(
    legend.position = "bottom", panel.grid.major.y = element_blank(),
    axis.text.y = element_text(
      colour = "transparent", face = "bold", size = 6.5, lineheight = 0.9, margin = margin(r = 3)
    ),
    strip.text = element_text(face = "bold")
  ) +
  guides(fill = guide_legend(nrow = 2))

save_panel(plot, "F01_Global", "e_protein_bars", 60, 62)
list(plot = plot, data = counts, inputs = inputs)
