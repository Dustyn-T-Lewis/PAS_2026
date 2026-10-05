#!/usr/bin/env Rscript
# F01 b. Proteins that differ between the sexes, in vehicle and in PAS, at FDR ≤ 0.05. Up and down
# are disjoint, so they share one bar.

source(here::here("R", "panels.R"))

inputs <- c(
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
counts <- readRDS(inputs[["contrasts"]])$results |>
  filter(universe == "full", contrast %in% c("sex_in_VEH", "sex_in_PAS"), adj.P.Val <= ALPHA) |>
  count(contrast, direction = if_else(logFC > 0, "Up", "Down")) |>
  mutate(
    contrast = factor(CONTRAST_LABELS[contrast], rev(CONTRAST_LABELS)),
    signed = if_else(direction == "Down", -n, n)
  )

panel <- ggplot(counts, aes(signed, contrast)) +
  contrast_bands() +
  geom_col(aes(fill = direction), width = 0.6) +
  geom_text(aes(label = n, hjust = if_else(direction == "Down", 1.2, -0.2)), size = 2) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = DIRECTION_COLOURS, name = NULL) +
  scale_x_continuous(labels = abs, expand = expansion(mult = 0.2)) +
  labs(
    x = "Proteins (lower | higher in males)", y = NULL, title = "Baseline sex differences",
    subtitle = "Male minus female, FDR ≤ 0.05"
  ) +
  theme_figure() +
  theme(legend.position = "none")
panel_data <- select(counts, contrast, direction, proteins = n)
panel_note <- "Proteins at FDR ≤ 0.05 per sex contrast and direction, full universe."
save_panel(panel, "F01", "main", "B_sex_proteins", 90, 50)
