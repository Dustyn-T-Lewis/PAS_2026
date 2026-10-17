# S5 a: features that put every PAS mouse above or below every vehicle mouse, per level and sex,
# against the same count under every relabelling of the sex's mice: median and 95th centile of the
# null, and the exact p (03_Classification).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(
  classification = output_path("02_Statistics", "03_Classification", "classification.rds")
)
LEVELS <- c(pathway = "Pathways", module = "Modules", complex = "Complexes", cluster = "Clusters")
perfect <- readRDS(inputs[["classification"]])$perfect |>
  mutate(
    level = factor(LEVELS[level], rev(LEVELS)),
    sex = factor(str_to_sentence(comparison), c("Female", "Male")),
    text = sprintf("p = %.2g", p)
  )

plot <- ggplot(perfect, aes(y = level)) +
  geom_segment(aes(x = 0, xend = null_q95), colour = "grey80", linewidth = 2.2) +
  geom_point(aes(x = null_median), shape = 124, size = 3, colour = "grey40") +
  geom_point(aes(x = observed, colour = sex), size = 1.8) +
  geom_text(aes(x = pmax(observed, null_q95), label = text),
    hjust = 0, nudge_x = 2.5, size = 1.7
  ) +
  facet_wrap(~sex) +
  scale_colour_manual(values = SEX_COLOURS, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.45))) +
  labs(
    x = "features separating perfectly", y = NULL, title = "Perfect separators against relabelling",
    subtitle = "Dot: observed; bar: 0 to the null's 95th centile; tick: null median"
  ) +
  theme_figure() +
  theme(panel.grid.major.y = element_blank())

save_panel(plot, "S5_Structure", "a_perfect_separators", 90, 55)
list(plot = plot, data = select(perfect, -text), inputs = inputs)
