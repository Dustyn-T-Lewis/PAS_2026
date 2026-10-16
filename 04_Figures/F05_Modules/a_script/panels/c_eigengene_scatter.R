# F05 c: the four strongest module-phenotype associations (pooled model, smallest p over both
# sexes): each mouse's eigengene against its trait value within the sex. The label gives Pearson r
# across these points and the p of the limma slope from 01_Fit, which also adjusts for the
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
    text = sprintf("Pearson r = %.2f\nslope p = %.2g", r, P.Value)
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
  geom_text(
    aes(-Inf, Inf, label = text), top,
    hjust = -0.08, vjust = 1.15, size = 1.7, lineheight = 0.9, inherit.aes = FALSE
  ) +
  facet_wrap(~cell, ncol = 2, scales = "free", labeller = label_wrap_gen(24)) +
  scale_colour_manual(values = GROUP_COLOURS, labels = GROUP_LABELS, name = NULL, drop = TRUE) +
  scale_x_continuous(breaks = scales::breaks_pretty(3)) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.3))) +
  labs(
    x = "trait value", y = "eigengene", title = "Strongest module–trait pairs",
    subtitle = "One point per mouse; slope p adjusted for SV1"
  ) +
  theme_figure() +
  theme(strip.text = element_text(size = 5), strip.clip = "off")

save_panel(plot, "F05_Modules", "c_eigengene_scatter", 70, 80)
list(
  plot = plot, data = left_join(points, select(top, cell, r, P.Value), by = "cell"),
  inputs = inputs
)
