#!/usr/bin/env Rscript
# F01 c. Proteins PAS moves in each sex and for the interaction at p ≤ 0.05 and Π ≤ 0.05. The two
# tiers overlap but are not nested, so they sit side by side.

source(here::here("R", "panels.R"))

TIERS <- c("p ≤ 0.05", "Π ≤ 0.05")

inputs <- c(
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
full <- readRDS(inputs[["contrasts"]])$results |>
  filter(universe == "full", contrast %in% PRIMARY) |>
  mutate(direction = if_else(logFC > 0, "Up", "Down"))
fdr <- count(filter(full, adj.P.Val <= ALPHA), contrast, direction)
counts <- full |>
  summarise(
    `p ≤ 0.05` = sum(P.Value <= ALPHA), `Π ≤ 0.05` = sum(pi_score <= ALPHA),
    .by = c(contrast, direction)
  ) |>
  pivot_longer(all_of(TIERS), names_to = "tier", values_to = "proteins") |>
  mutate(
    contrast = factor(CONTRAST_LABELS[contrast], rev(CONTRAST_LABELS[PRIMARY])),
    tier = factor(tier, rev(TIERS)),
    signed = if_else(direction == "Down", -proteins, proteins)
  )
chance <- round(n_distinct(full$protein) * ALPHA / 2)
fdr_note <- fdr |>
  summarise(n = sum(n), .by = contrast) |>
  mutate(text = str_glue("{n} in {str_remove(CONTRAST_LABELS[contrast], 'PAS vs VEH, ')}"))

panel <- ggplot(counts, aes(signed, contrast, group = tier)) +
  contrast_bands() +
  geom_vline(xintercept = c(-chance, chance), linetype = "dashed", linewidth = 0.3) +
  geom_col(
    aes(fill = direction, alpha = tier),
    position = position_dodge(width = 0.75), width = 0.7
  ) +
  geom_text(
    aes(label = proteins, hjust = if_else(direction == "Down", 1.15, -0.15)),
    position = position_dodge(width = 0.75), size = 1.6
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = DIRECTION_COLOURS, guide = "none") +
  scale_alpha_manual(values = c("p ≤ 0.05" = 0.4, "Π ≤ 0.05" = 1), breaks = TIERS, name = NULL) +
  scale_x_continuous(labels = abs, expand = expansion(mult = 0.2)) +
  labs(
    x = "Proteins (down | up)", y = NULL, title = "Proteins PAS moves",
    subtitle = str_glue(
      "Π = p^|log2 FC| ranks, it controls no error rate; dashed: {chance} by chance;\n",
      "FDR ≤ 0.05: {paste(fdr_note$text, collapse = ', ')}; none in males or for the interaction"
    )
  ) +
  theme_figure() +
  theme(legend.position = "bottom", legend.margin = margin(t = -4))
panel_data <- select(counts, contrast, direction, tier, proteins)
panel_note <- "Proteins at p ≤ 0.05 and Π ≤ 0.05 per PAS contrast and direction, full universe."
save_panel(panel, "F01", "main", "C_pas_proteins", 90, 65)
