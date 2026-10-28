# S1 a: precursors and protein groups entering filtering (past the three q-value filters) and left
# after each filter, from 01_Filtering.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

# 01_Filtering's step names, as the panel prints them.
STEP_LABELS <- c(
  "q-value filters passed" = "passed q-value filters",
  "cRAP groups removed" = "cRAP-only groups removed",
  "shared peptides removed" = "shared peptides removed",
  "compound groups removed" = "compound groups removed",
  "contaminants removed" = "keratin, Ig, haemoglobin removed"
)

inputs <- c(filtering = output_path("01_Processing", "01_Filtering", "precursors_filtered.rds"))
steps <- readRDS(inputs[["filtering"]])$filter_log |>
  mutate(step = factor(step, rev(step))) |>
  pivot_longer(c(precursors, protein_groups), names_to = "unit", values_to = "n") |>
  mutate(unit = factor(
    if_else(unit == "precursors", "Precursors", "Protein groups"),
    c("Precursors", "Protein groups")
  ))

plot <- ggplot(steps, aes(n, step)) +
  geom_col(fill = "grey70", width = 0.65) +
  geom_text(aes(label = format(n, big.mark = ",")), hjust = 1.08, size = 1.7, colour = "white") +
  facet_wrap(~unit, scales = "free_x") +
  scale_x_continuous(
    labels = scales::label_comma(), breaks = scales::breaks_pretty(2),
    expand = expansion(mult = c(0, 0.04))
  ) +
  scale_y_discrete(labels = STEP_LABELS) +
  labs(
    x = "count", y = NULL, title = "Filtering", subtitle = "Counts entering, then after each filter"
  ) +
  theme_figure() +
  theme(panel.grid.major.y = element_blank(), panel.spacing.x = unit(4, "mm"))

save_panel(plot, "S1_Processing", "a_filtering", 70, 55)
list(plot = plot, data = steps, inputs = inputs)
