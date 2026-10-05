#!/usr/bin/env Rscript
# S2 c. Distribution of log2 fold changes per contrast, full universe.

source(here::here("R", "panels.R"))

inputs <- c(
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
full <- readRDS(inputs[["contrasts"]])$results |>
  filter(universe == "full") |>
  mutate(contrast = factor(CONTRAST_LABELS[contrast], rev(CONTRAST_LABELS)))
medians <- summarise(full, median_abs = median(abs(logFC)), .by = contrast)
panel <- ggplot(full, aes(logFC, contrast, fill = contrast)) +
  ggridges::geom_density_ridges(
    colour = "grey30", linewidth = 0.25, scale = 1.2, rel_min_height = 0.005, alpha = 0.8
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_text(
    data = medians, aes(x = Inf, label = sprintf("%.2f", median_abs)),
    hjust = 1.1, vjust = -0.6, size = 1.9
  ) +
  scale_fill_manual(values = CONTRAST_COLOURS, guide = "none") +
  coord_cartesian(xlim = c(-1.5, 1.5)) +
  labs(
    x = "log2 fold change", y = NULL, title = "Effect sizes",
    subtitle = "Median |log2 FC| at right"
  ) +
  theme_figure()
panel_data <- medians
panel_note <- "Median absolute log2 fold change per contrast, full universe."
save_panel(panel, "F01", "supp", "S2_C_effect_sizes", 80, 65)
