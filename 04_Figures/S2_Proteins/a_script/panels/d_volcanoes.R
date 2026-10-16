# S2 d: every protein in the full proteome per contrast, log2 FC against -log10 p, Π <= 0.05
# coloured by direction and FDR <= 0.05 outlined (01_Fit).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(fit = output_path("02_Statistics", "01_Fit", "fit.rds"))
proteins <- readRDS(inputs[["fit"]])$protein_contrasts |>
  filter(universe == "full", contrast %in% COMPARISON_CONTRASTS) |>
  mutate(
    contrast = factor(CONTRAST_SHORT[contrast], CONTRAST_SHORT),
    status = case_when(
      pi_score > ALPHA ~ "Π > 0.05", logFC > 0 ~ "Up, Π ≤ 0.05",
      .default = "Down, Π ≤ 0.05"
    ),
    status = factor(status, c("Down, Π ≤ 0.05", "Up, Π ≤ 0.05", "Π > 0.05"))
  ) |>
  arrange(desc(status))

plot <- ggplot(proteins, aes(logFC, -log10(P.Value))) +
  geom_point(aes(colour = status), size = 0.5, alpha = 0.8) +
  geom_point(data = filter(proteins, adj.P.Val <= ALPHA), shape = 21, size = 1.6, stroke = 0.4) +
  geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", linewidth = 0.25) +
  facet_wrap(~contrast, nrow = 1) +
  scale_colour_manual(
    values = c(
      "Down, Π ≤ 0.05" = DIRECTION_COLOURS[["Down"]], "Up, Π ≤ 0.05" = DIRECTION_COLOURS[["Up"]],
      "Π > 0.05" = "grey80"
    ),
    name = NULL, guide = guide_legend(override.aes = list(size = 1.5))
  ) +
  labs(
    x = "log2 fold change", y = "-log10 p", title = "Every protein per contrast",
    subtitle = "Full proteome; dashed: p = 0.05; ringed: FDR ≤ 0.05"
  ) +
  theme_figure()

save_panel(plot, "S2_Proteins", "d_volcanoes", 180, 62)
list(
  plot = plot,
  data = select(proteins, contrast, feature, gene, logFC, P.Value, adj.P.Val, pi_score, status),
  inputs = inputs
)
