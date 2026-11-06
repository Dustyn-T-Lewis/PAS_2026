# S1 f: what went into each pellet beyond mitochondria. Per mouse, the mean scaled log2 abundance
# of marker proteins for sarcoplasmic and endoplasmic reticulum, myofibril, nerve, red cell and
# vessel, and haemoglobin's share of signal, with the PAS p within each sex (05_Composition).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(composition = output_path("02_Statistics", "05_Composition", "composition.rds"))
comp <- readRDS(inputs[["composition"]])
MARKERS <- c(comp$groups$markers, comp$groups$haemoglobin)
shown <- filter(comp$composition, group %in% MARKERS)
# Each column's unit as a strip above its name.
UNITS <- c("marker score", "% of signal")
facets <- \(data) {
  mutate(data,
    unit = factor(if_else(group %in% comp$groups$markers, UNITS[1], UNITS[2]), UNITS),
    marker = factor(group, MARKERS, c(comp$groups$markers, "Haemoglobin"))
  )
}
scores <- shown |>
  facets() |>
  mutate(group = factor(group_label, FIGURE_GROUPS))
labels <- comp$tests |>
  filter(group %in% MARKERS) |>
  summarise(
    text = paste(paste0(names(SEX_NAMES)[match(sex, SEX_NAMES)], ": ", format_p(p)), collapse = "\n"),
    .by = group
  ) |>
  facets()

plot <- ggplot(scores, aes(group, score)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.25, colour = "grey30") +
  group_points(position_jitter(width = 0.12, seed = SEED)) +
  geom_text(
    aes(-Inf, Inf, label = text), labels,
    hjust = -0.08, vjust = 1.2, size = 1.6, lineheight = 0.9, inherit.aes = FALSE
  ) +
  ggh4x::facet_nested_wrap(~ unit + marker,
    nrow = 1, scales = "free_y", labeller = labeller(marker = label_wrap_gen(14))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(
    x = NULL, y = NULL, title = "Pellet composition",
    subtitle = paste(
      "Marker score: mean scaled log2 abundance of a group's marker proteins;",
      "haemoglobin: % of precursor signal"
    )
  ) +
  theme_figure() +
  theme(axis.text.x = element_blank(), strip.text = element_text(size = 5.5))

save_figure(plot, "F01_Global", "S1f_composition", 180, 60)
list(plot = plot, data = rename(shown, marker = group, group = group_label))
