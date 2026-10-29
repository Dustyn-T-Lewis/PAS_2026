# S1 b: precursors and proteins identified per run, by group, from 01_Filtering.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(filtering = output_path("01_Processing", "01_Filtering", "precursors_filtered.rds"))
runs <- readRDS(inputs[["filtering"]])$run_qc |>
  select(run, sample_id, group, precursors, proteins) |>
  pivot_longer(c(precursors, proteins), names_to = "unit", values_to = "n") |>
  mutate(
    group = factor(group, FIGURE_GROUPS),
    unit = factor(str_to_sentence(unit), c("Precursors", "Proteins"))
  )

plot <- ggplot(runs, aes(group, n)) +
  geom_point(aes(fill = group, shape = str_sub(group, 1, 3)),
    position = position_jitter(width = 0.15, seed = SEED), size = 1.5, colour = "grey15",
    stroke = 0.25
  ) +
  facet_wrap(~unit, scales = "free_y") +
  group_point_scales() +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x = NULL, y = "identified per run", title = "Depth per run",
    subtitle = "One point per mouse"
  ) +
  theme_figure() +
  theme(axis.text.x = element_blank())

save_figure(plot, "F01_Global", "S1b_run_depth", 55, 55)
list(plot = plot, data = runs)
