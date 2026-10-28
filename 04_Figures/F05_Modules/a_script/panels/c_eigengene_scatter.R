# F05 c: the four strongest module-phenotype associations (pooled model, smallest p over both
# sexes): each mouse's eigengene against its trait value within the sex. The facet title gives
# Pearson r across these points and the p of the limma slope from 01_Fit, which also adjusts for the
# surrogate variable; the two can differ.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(fit = output_path("02_Statistics", "01_Fit", "fit.rds"))
fit <- readRDS(inputs[["fit"]])

top <- fit$score_associations |>
  filter(level == "module", model == "pooled") |>
  slice_min(P.Value, n = 4, with_ties = FALSE) |>
  mutate(
    sex_code = if_else(sex == "female", "F", "M"),
    cell = sprintf("%s ~ %s, %s", feature, PHENOTYPE_LABELS[phenotype], sex),
    cell = factor(cell, cell),
    strip = sprintf("%s\nPearson r = %.2f, slope p = %.2g", str_wrap(cell, 24), r, P.Value)
  )
eigengenes <- as_tibble(t(fit$scores$module), rownames = "sample_id")
mice <- fit$targets |>
  select(sample_id, mouse, sex, group) |>
  left_join(select(fit$phenotypes, mouse, all_of(names(PHENOTYPE_LABELS))), by = "mouse") |>
  left_join(eigengenes, by = "sample_id")
points <- top |>
  select(feature, phenotype, sex_code, cell) |>
  pmap(\(feature, phenotype, sex_code, cell) {
    mice |>
      filter(sex == sex_code) |>
      transmute(cell, sample_id, group, trait = .data[[phenotype]], eigengene = .data[[feature]])
  }) |>
  list_rbind() |>
  mutate(group = factor(group, FIGURE_GROUPS))

plot <- ggplot(points, aes(trait, eigengene)) +
  geom_point(aes(colour = group), size = 1.3) +
  facet_wrap(
    ~cell,
    ncol = 2, scales = "free", labeller = as_labeller(deframe(select(top, cell, strip)))
  ) +
  scale_colour_manual(values = GROUP_COLOURS, labels = GROUP_LABELS, name = NULL, drop = TRUE) +
  scale_x_continuous(breaks = scales::breaks_pretty(3)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  labs(
    x = "trait value", y = "eigengene", title = "Strongest module–trait pairs",
    subtitle = "One point per mouse; slope p adjusted for SV1"
  ) +
  theme_figure() +
  theme(strip.text = element_text(size = 5), strip.clip = "off")

save_panel(plot, "F05_Modules", "c_eigengene_scatter", 110, 68)
list(
  plot = plot, data = left_join(points, select(top, cell, r, P.Value), by = "cell"),
  inputs = inputs
)
