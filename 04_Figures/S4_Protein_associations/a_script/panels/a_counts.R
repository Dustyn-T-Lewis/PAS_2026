# S4 a: proteins associated with each phenotype at p <= 0.05 within each sex, under each
# association model (01_Fit, full proteome), split by the sign of the limma t: positive to the
# right, negative to the left.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(fit = output_path("02_Statistics", "01_Fit", "fit.rds"))
DIRECTION_KEYS <- c(Up = "▲ positive", Down = "▼ negative")
counts <- readRDS(inputs[["fit"]])$protein_associations |>
  filter(universe == "full") |>
  mutate(direction = if_else(t > 0, "Up", "Down")) |>
  summarise(
    p_le_0.05 = sum(P.Value <= ALPHA), fdr_le_0.05 = sum(adj.P.Val <= ALPHA),
    .by = c(phenotype, sex, model, direction)
  ) |>
  complete(phenotype, sex, model, direction, fill = list(p_le_0.05 = 0L, fdr_le_0.05 = 0L)) |>
  mutate(
    trait = factor(PHENOTYPE_LABELS[phenotype], rev(PHENOTYPE_LABELS)),
    sex = factor(str_to_sentence(sex), c("Female", "Male")),
    model = factor(MODEL_LABELS[model], MODEL_LABELS),
    direction = factor(direction, names(DIRECTION_KEYS)),
    x = if_else(direction == "Up", 1L, -1L) * p_le_0.05
  ) |>
  arrange(model, sex, trait, direction)
pad <- 0.03 * max(counts$p_le_0.05)

plot <- ggplot(counts, aes(x, trait, fill = direction)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.25) +
  geom_text(
    aes(x + sign(x) * pad,
      label = p_le_0.05, colour = direction,
      hjust = if_else(direction == "Up", 0, 1)
    ),
    size = 1.7
  ) +
  facet_grid(sex ~ model) +
  scale_fill_manual(values = DIRECTION_COLOURS, labels = DIRECTION_KEYS, name = NULL) +
  scale_colour_manual(values = DIRECTION_COLOURS, guide = "none") +
  scale_x_continuous(labels = abs, expand = expansion(mult = 0.16)) +
  labs(
    x = "proteins at p ≤ 0.05", y = NULL, title = "Protein–phenotype associations",
    subtitle = "Full proteome, within sex; direction from the sign of the limma t"
  ) +
  theme_figure() +
  theme(panel.grid.major.y = element_blank(), panel.spacing.x = unit(3, "mm"))

save_panel(plot, "S4_Protein_associations", "a_counts", 180, 62)
list(
  plot = plot,
  data = select(counts, phenotype, sex, model, direction, p_le_0.05, fdr_le_0.05),
  inputs = inputs
)
