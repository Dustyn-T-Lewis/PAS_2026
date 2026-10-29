# S5 b: every protein at FDR ≤ 0.10 in the baseline, its FDR under each variant. Tile colour bins
# the FDR; the text is the log2 FC, so a sign change would show.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(sensitivity = output_path("02_Statistics", "07_Sensitivity", "sensitivity.rds"))
sensitivity <- readRDS(inputs[["sensitivity"]])
named <- sensitivity$named
FDR_BINS <- c("≤ 0.05", "≤ 0.10", "≤ 0.25", "> 0.25")

order <- named |>
  filter(variant == "baseline") |>
  arrange(contrast, adj.P.Val)
cells <- named |>
  mutate(
    variant = factor(variant, levels(sensitivity$counts$variant)),
    contrast = factor(contrast, COMPARISON_CONTRASTS),
    row = factor(paste(contrast, feature), rev(paste(order$contrast, order$feature))),
    fdr = cut(adj.P.Val, c(0, 0.05, 0.10, 0.25, 1), FDR_BINS, include.lowest = TRUE)
  )
genes <- set_names(order$gene, paste(order$contrast, order$feature))

plot <- ggplot(cells, aes(variant, row, fill = fdr)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", logFC), colour = fdr == FDR_BINS[1]), size = 1.7) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  ggh4x::facet_grid2(contrast ~ .,
    scales = "free_y", space = "free_y", labeller = as_labeller(CONTRAST_SHORT),
    strip = contrast_strips(y = levels(droplevels(cells$contrast)))
  ) +
  scale_y_discrete(labels = genes) +
  scale_x_discrete(labels = scales::label_wrap(12)) +
  scale_fill_manual(
    values = set_names(c("#4D4D4D", "#999999", "#D9D9D9", "#F5F5F5"), FDR_BINS),
    name = "FDR", drop = FALSE
  ) +
  labs(
    x = NULL, y = NULL, title = "Baseline proteins under each variant",
    subtitle = "Proteins at FDR ≤ 0.10 in the baseline; text: log2 FC"
  ) +
  theme_figure() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.y = element_text(face = "italic")
  ) +
  guides(fill = guide_legend(nrow = 1))

save_figure(plot, "F02_Proteins", "S5b_named", 180, 100)
list(plot = plot, data = named)
