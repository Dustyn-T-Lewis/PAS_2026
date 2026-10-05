#!/usr/bin/env Rscript
# F01 a. PCA of the twenty mice: group ellipses, PERMANOVA, and the MitoCarta share as a fitted
# vector (05_Summary/05_Global).

source(here::here("R", "panels.R"))

inputs <- c(global = here("05_Summary", "05_Global", "c_data", "global.rds"))
global <- readRDS(inputs[["global"]])
scores <- mutate(global$scores, group = factor(group, FIGURE_GROUPS))
stat <- \(t) with(filter(global$permanova, term == t), sprintf("R² %.2f, p = %.3g", r2, p))
reach <- 0.8 * min(diff(range(scores$PC1)), diff(range(scores$PC2))) / 2
arrow <- mutate(global$mito_vector, x = PC1 * sqrt(r2) * reach, y = PC2 * sqrt(r2) * reach)

panel <- ggplot(scores, aes(PC1, PC2)) +
  stat_ellipse(aes(fill = group), geom = "polygon", alpha = 0.18, colour = NA) +
  stat_ellipse(aes(colour = group), linewidth = 0.3) +
  geom_point(aes(colour = group, shape = sex), size = 2) +
  geom_segment(
    data = arrow, aes(x = 0, y = 0, xend = x, yend = y), inherit.aes = FALSE,
    arrow = grid::arrow(length = unit(1.5, "mm")), linewidth = 0.5
  ) +
  geom_label(
    data = arrow, aes(x = x, y = y, label = sprintf("MitoCarta share\nR² %.2f, p = %.3g", r2, p)),
    inherit.aes = FALSE, hjust = 0.5, vjust = 1.3, size = 1.9, lineheight = 0.9,
    fill = alpha("white", 0.85), linewidth = 0, label.padding = unit(0.8, "mm")
  ) +
  scale_colour_manual(values = GROUP_COLOURS, name = NULL) +
  scale_fill_manual(values = GROUP_COLOURS, guide = "none") +
  scale_shape_manual(values = SEX_SHAPES, guide = "none") +
  labs(
    x = sprintf("PC1 (%.1f%%)", 100 * global$variance[1]),
    y = sprintf("PC2 (%.1f%%)", 100 * global$variance[2]),
    title = "Proteome ordination",
    subtitle = str_glue(
      "PERMANOVA: sex {stat('sex')}; treatment {stat('treatment')};\n",
      "interaction {stat('sex:treatment')}; MitoCarta share fitted with envfit"
    )
  ) +
  theme_figure() +
  theme(legend.position = "inside", legend.position.inside = c(0.12, 0.85))
panel_data <- scores
panel_note <- "PCA scores per mouse (05_Summary/05_Global)."
save_panel(panel, "F01", "main", "A_pca", 110, 95)
