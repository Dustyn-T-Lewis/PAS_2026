# S1 b: precursors and proteins identified per run, by group, from 01_Filtering.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(filtering = output_path("01_Processing", "01_Filtering", "01_filtering.xlsx"))
runs <- readxl::read_excel(inputs[["filtering"]], "run_qc") |>
  select(run, sample_id, group, precursors, proteins) |>
  pivot_longer(c(precursors, proteins), names_to = "unit", values_to = "n") |>
  mutate(
    group = factor(group, FIGURE_GROUPS),
    unit = factor(str_to_sentence(unit), c("Precursors", "Proteins"))
  )

plot <- ggplot(runs, aes(group, n, colour = group)) +
  geom_point(position = position_jitter(width = 0.15, seed = SEED), size = 1.2) +
  facet_wrap(~unit, scales = "free_y") +
  scale_colour_manual(values = GROUP_COLOURS, labels = GROUP_LABELS, guide = "none") +
  scale_x_discrete(labels = GROUP_SHORT) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x = NULL, y = "identified per run", title = "Depth per run",
    subtitle = "One point per mouse; groups coloured as in d"
  ) +
  theme_figure() +
  theme(axis.text.x = element_blank(), strip.text = element_text(face = "bold"))

save_panel(plot, "S1_Processing", "b_run_depth", 55, 55)
list(plot = plot, data = runs, inputs = inputs)
