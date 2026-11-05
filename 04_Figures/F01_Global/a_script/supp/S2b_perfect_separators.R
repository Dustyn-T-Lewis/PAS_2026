# S2 b: features that put every PAS mouse above or below every vehicle mouse, per level and sex,
# against the same count under every relabelling of the sex's mice: median and 95th centile of the
# null, and the exact p (03_Classification).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(
  classification = output_path("02_Statistics", "03_Classification", "classification.rds")
)
LEVELS <- c(pathway = "Pathways", module = "Modules", complex = "Complexes", cluster = "Clusters")
perfect <- readRDS(inputs[["classification"]])$perfect |>
  mutate(
    level = factor(LEVELS[level], LEVELS),
    sex = factor(str_to_sentence(comparison), c("Male", "Female"))
  )

# Each level gets its own x range, shared by the sexes so they compare; the p text sits a fixed
# share of that range right of the longer of bar and dot, clear of the dot, and an invisible point
# widens the range by ROOM to fit it. Breaks stop at the data so no tick sits in the text's room.
ROOM <- 1.6
labelled <- perfect |>
  mutate(
    end = pmax(observed, null_q95),
    range = max(end, 1),
    text_x = end + 0.1 * range,
    room = ROOM * range,
    text = sprintf("p = %.2g", p),
    .by = level
  )
integer_breaks <- \(limits) {
  top <- limits[2] / ROOM
  breaks <- scales::breaks_pretty(n = 3)(c(0, top))
  breaks[breaks == round(breaks) & breaks <= top]
}

plot <- ggplot(labelled, aes(y = sex)) +
  geom_segment(aes(x = 0, xend = null_q95), colour = "grey80", linewidth = 2.2) +
  geom_point(aes(x = null_median), shape = 124, size = 3, colour = "grey40") +
  geom_point(aes(x = observed, colour = sex), size = 1.8) +
  geom_text(aes(x = text_x, label = text), hjust = 0, size = 1.7) +
  geom_blank(aes(x = room)) +
  facet_wrap(~level, ncol = 2, scales = "free_x") +
  scale_colour_manual(values = SEX_COLOURS, guide = "none") +
  scale_x_continuous(breaks = integer_breaks, expand = expansion(mult = c(0.04, 0))) +
  labs(
    x = "features separating perfectly", y = NULL, title = "Perfect separators against relabelling",
    subtitle = "Dot: observed; bar: 0 to the null's 95th centile; tick: null median"
  ) +
  theme_figure() +
  theme(panel.grid.major.y = element_blank(), panel.spacing.x = unit(4, "mm"))

save_figure(plot, "F01_Global", "S2b_perfect_separators", 90, 55)
list(plot = plot, data = select(perfect, -sex))
