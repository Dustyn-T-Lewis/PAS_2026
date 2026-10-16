# S1 f: what went into each pellet beyond mitochondria. Per mouse, the mean scaled log2 abundance
# of marker proteins for sarcoplasmic and endoplasmic reticulum, myofibril, nerve, red cell and
# vessel, and haemoglobin's share of signal, with the PAS p within each sex (05_Composition).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(composition = output_path("02_Statistics", "05_Composition", "composition.rds"))
comp <- readRDS(inputs[["composition"]])
MARKERS <- c(
  "Sarcoplasmic reticulum", "Endoplasmic reticulum", "Myofibril", "Nerve", "Red cell", "Vessel",
  "Haemoglobin (% of signal)"
)
scores <- comp$composition |>
  filter(group %in% MARKERS) |>
  mutate(marker = factor(group, MARKERS), group = factor(group_label, FIGURE_GROUPS))
labels <- comp$tests |>
  filter(group %in% MARKERS) |>
  summarise(
    text = paste(sprintf("%s: p = %.2g", if_else(sex == "female", "F", "M"), p), collapse = "\n"),
    .by = group
  ) |>
  mutate(marker = factor(group, MARKERS))

plot <- ggplot(scores, aes(group, score, colour = group)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.25, colour = "grey30") +
  geom_point(position = position_jitter(width = 0.12, seed = SEED), size = 0.9) +
  geom_text(
    aes(-Inf, Inf, label = text), labels,
    hjust = -0.08, vjust = 1.2, size = 1.6, lineheight = 0.9, inherit.aes = FALSE
  ) +
  facet_wrap(~marker, nrow = 1, scales = "free_y", labeller = label_wrap_gen(14)) +
  scale_colour_manual(values = GROUP_COLOURS, labels = GROUP_LABELS, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(
    x = NULL, y = "score", title = "Pellet composition",
    subtitle = "Mean scaled log2 abundance of marker proteins; haemoglobin: % of raw signal"
  ) +
  theme_figure() +
  theme(axis.text.x = element_blank(), strip.text = element_text(face = "bold", size = 5.5)) +
  guides(colour = guide_legend(nrow = 1))

save_panel(plot, "S1_Processing", "f_composition", 180, 60)
list(plot = plot, data = comp$composition, inputs = inputs)
