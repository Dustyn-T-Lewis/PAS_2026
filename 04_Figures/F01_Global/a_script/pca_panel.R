# The PCA panels F01 b and c share, sourced by both.

# PCA scores with a shaded 95% t-ellipse per group (fill: group; shape: treatment), the PERMANOVA
# terms in the top-right corner (bold at p <= 0.05), and arrows (x, y, label) from the origin, each
# in its own colour with a key. An arrow of length 1 on an axis reaches 80% of the largest score.
pca_panel <- function(scores, variance, title, subtitle, permanova, arrows, arrow_colours) {
  axis <- \(i) sprintf("PC%d (%.1f%%)", i, 100 * variance[i])
  scores <- mutate(scores, group = factor(group, FIGURE_GROUPS))
  reach <- 0.8 * max(abs(c(scores$PC1, scores$PC2)))
  arrows <- mutate(arrows,
    x = x * reach, y = y * reach, label = factor(label, names(arrow_colours))
  )
  terms <- c(sex = "Sex", treatment = "Treatment", `sex:treatment` = "Treatment × sex")
  permanova <- permanova |>
    filter(term %in% names(terms)) |>
    transmute(
      text = sprintf("%s  R² %.2f, p %.2g", terms[term], r2, p),
      face = if_else(p <= ALPHA, "bold", "plain"),
      vjust = 1.5 + 1.35 * (row_number() - 1)
    )
  ggplot(scores, aes(PC1, PC2)) +
    stat_ellipse(
      aes(fill = group, colour = after_scale(fill)),
      geom = "polygon", type = "t", level = 0.95, linewidth = 0.3, alpha = 0.12,
      show.legend = FALSE
    ) +
    geom_segment(
      aes(0, 0, xend = x, yend = y, colour = label), arrows,
      arrow = grid::arrow(length = grid::unit(1.1, "mm"), type = "closed"), linewidth = 0.45
    ) +
    group_points() +
    geom_text(
      aes(Inf, Inf, label = text, fontface = face, vjust = vjust), permanova,
      hjust = 1.04, size = 1.6, colour = "grey15", inherit.aes = FALSE
    ) +
    scale_colour_manual(
      values = arrow_colours, labels = scales::label_wrap(20), name = NULL,
      guide = guide_legend(ncol = 2, order = 2, override.aes = list(linewidth = 0.5))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    labs(x = axis(1), y = axis(2), title = title, subtitle = subtitle) +
    theme_figure() +
    theme(legend.box = "vertical", legend.key.width = unit(4, "mm"))
}
