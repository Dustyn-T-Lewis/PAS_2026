#!/usr/bin/env Rscript
# S2 e. Female against male PAS log2 fold change per protein, by universe.

source(here::here("R", "panels.R"))

UNIVERSE_LABELS <- c(full = "Full proteome", mito = "Mitochondrial universe")
inputs <- c(
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
effects <- readRDS(inputs[["contrasts"]])$results |>
  filter(contrast %in% c("PAS_vs_VEH_female", "PAS_vs_VEH_male")) |>
  select(universe, protein, contrast, logFC) |>
  pivot_wider(names_from = contrast, values_from = logFC) |>
  mutate(universe = factor(UNIVERSE_LABELS[universe], UNIVERSE_LABELS))
correlation <- effects |>
  summarise(
    rho = cor(PAS_vs_VEH_female, PAS_vs_VEH_male, method = "spearman"), proteins = n(),
    .by = universe
  ) |>
  mutate(label = sprintf("ρ = %.2f\n%d proteins", rho, proteins))
panel <- ggplot(effects, aes(PAS_vs_VEH_male, PAS_vs_VEH_female)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey60") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey60") +
  geom_point(size = 0.35, alpha = 0.5, colour = "grey30") +
  geom_text(
    data = correlation, aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.2, size = 1.9
  ) +
  facet_wrap(~universe) +
  labs(
    x = "log2 FC, PAS vs VEH, males", y = "log2 FC, females",
    title = "Female against male PAS effect",
    subtitle = "Spearman ρ; descriptive, the interaction contrast is the test"
  ) +
  theme_figure() +
  theme(strip.text = element_text(size = 6))
panel_data <- select(correlation, -label)
panel_note <- "Spearman correlation of female and male PAS log2 FC per universe."
save_panel(panel, "F01", "supp", "S2_E_female_male", 110, 65)
