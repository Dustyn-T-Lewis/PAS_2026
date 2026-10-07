# Figure style and the plots more than one script draws, for figures and substage reports alike.
# Sourced by R/helpers.R.

library(ggplot2)
library(patchwork)

PORTRAIT_MM <- 180
LANDSCAPE_MM <- 297
TREATMENT_COLOURS <- c(VEH = "#56B4E9", PAS = "#D55E00")
SEX_COLOURS <- c(Female = "#D55E00", Male = "#0072B2")
SEX_SHAPES <- c(F = 16, M = 17)
FIGURE_GROUPS <- c("VEH_F", "PAS_F", "VEH_M", "PAS_M")
GROUP_COLOURS <- c(VEH_F = "#FDB863", PAS_F = "#D55E00", VEH_M = "#92C5DE", PAS_M = "#0072B2")
DIRECTION_COLOURS <- c(Down = "#2166AC", Up = "#B2182B")
PHENOTYPE_LABELS <- c(
  csa_mean = "CSA", capillary_density = "Capillary density",
  satellite_cells_per_fiber = "Satellite cells / fibre",
  conductance_sol_pyr = "Soleus pyruvate", conductance_pla_pyr = "Plantaris pyruvate",
  conductance_sol_oct = "Soleus octanoyl"
)
CONTRAST_LABELS <- c(
  PAS_vs_VEH_female = "PAS vs VEH, female", PAS_vs_VEH_male = "PAS vs VEH, male",
  treatment_by_sex = "Treatment × sex", sex_in_VEH = "Sex in VEH", sex_in_PAS = "Sex in PAS"
)
CONTRAST_COLOURS <- set_names(
  c("#D55E00", "#0072B2", "#7B3294", "grey45", "grey70"), CONTRAST_LABELS
)
DATABASE_COLOURS <- c("GO:BP" = "#999999", Reactome = "#56B4E9", MitoCarta = "#009E73")

# Framed panels with a light grid, bold title and italic subtitle, after the YvO and HRvLR figures.
theme_figure <- function() {
  theme_bw(base_size = 7, base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = 7.5, face = "bold"),
      plot.subtitle = element_text(size = 5.5, face = "italic", colour = "grey35"),
      axis.text = element_text(size = 6),
      plot.tag = element_text(size = 9, face = "bold"),
      panel.border = element_rect(colour = "grey35", linewidth = 0.4),
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 7, face = "bold"),
      legend.key.size = unit(3, "mm"),
      legend.background = element_blank()
    )
}

# A faint band in each contrast's colour behind its row (or column), as in the YvO overview
# figures. Draw it first; it reads the contrast (and any facet column) from the panel's own data.
contrast_bands <- function(alpha = 0.12, axis = c("y", "x")) {
  axis <- match.arg(axis)
  band_data <- \(d) {
    distinct(d, across(any_of(c("contrast", "family")))) |>
      mutate(band = CONTRAST_COLOURS[as.character(contrast)])
  }
  if (axis == "y") {
    geom_tile(
      data = band_data, aes(x = 0, y = contrast, fill = I(band)),
      width = Inf, height = 1, alpha = alpha, inherit.aes = FALSE
    )
  } else {
    geom_tile(
      data = band_data, aes(x = contrast, y = 0, fill = I(band)),
      width = 1, height = Inf, alpha = alpha, inherit.aes = FALSE
    )
  }
}

save_figure <- function(figure, dir, name, height_mm, width_mm = PORTRAIT_MM) {
  save <- \(ext, ...) {
    ggsave(file.path(dir, paste0(name, ext)), figure,
      width = width_mm, height = height_mm, units = "mm", ...
    )
  }
  save(".pdf", device = grDevices::cairo_pdf)
  save(".png", dpi = 600)
}

# MDS of the mice from limpa, which weights distances by quantification standard errors.
mds_plot <- function(proteins, targets, title) {
  mds <- local({
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    limpa::plotMDSUsingSEs(proteins)
  })
  coords <- tibble(sample_id = colnames(proteins$E), dim1 = mds$x, dim2 = mds$y) |>
    left_join(select(targets, sample_id, sex, group), by = "sample_id") |>
    mutate(group = factor(group, FIGURE_GROUPS))
  ggplot(coords, aes(dim1, dim2, colour = group, shape = sex)) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = GROUP_COLOURS, name = NULL) +
    scale_shape_manual(values = SEX_SHAPES, guide = "none") +
    labs(
      x = sprintf("Dimension 1 (%.0f%%)", 100 * mds$var.explained[1]),
      y = sprintf("Dimension 2 (%.0f%%)", 100 * mds$var.explained[2]),
      title = title
    ) +
    theme_figure()
}

# One value or several per mouse, grouped and coloured by group, median marked.
per_mouse_plot <- function(data, targets, measures, title, subtitle, y) {
  data |>
    select(sample_id, all_of(names(measures))) |>
    pivot_longer(-sample_id, names_to = "measure", values_to = "value") |>
    mutate(measure = factor(measures[measure], measures)) |>
    left_join(select(targets, sample_id, sex, group), by = "sample_id") |>
    mutate(group = factor(group, FIGURE_GROUPS)) |>
    ggplot(aes(group, value, colour = group, shape = sex)) +
    stat_summary(
      fun = median, geom = "crossbar", width = 0.5, linewidth = 0.25,
      colour = "grey40"
    ) +
    geom_point(position = position_jitter(width = 0.12, seed = SEED), size = 1.4) +
    facet_wrap(~measure, scales = "free_y") +
    scale_colour_manual(values = GROUP_COLOURS) +
    scale_shape_manual(values = SEX_SHAPES) +
    labs(x = NULL, y = y, title = title, subtitle = subtitle) +
    theme_figure() +
    theme(
      legend.position = "none", axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
      strip.text = element_text(size = 6)
    )
}

# A substage report: one landscape page per named ggplot, cairo so Π and ≤ render. Returns the
# rows of the workbook's plot_index sheet.
write_report <- function(file, pages, width = 7, height = 5) {
  cairo_pdf(file, width = width, height = height, onefile = TRUE)
  on.exit(dev.off())
  walk(pages, print)
  tibble(file = basename(file), page = seq_along(pages), shows = names(pages))
}

# Rows z-scored, ordered by clustering on 1 - Pearson r (display only), mice in group order.
group_heatmap <- function(m, targets, title, subtitle = NULL, row_labels = nrow(m) <= 60) {
  z <- t(scale(t(m)))
  rows <- rownames(z)[cluster_order(z)]
  limit <- unname(quantile(abs(z), 0.99, na.rm = TRUE))
  as_tibble(z, rownames = "feature") |>
    pivot_longer(-feature, names_to = "sample_id", values_to = "z") |>
    left_join(select(targets, sample_id, group), by = "sample_id") |>
    mutate(
      feature = factor(feature, rows), group = factor(group, FIGURE_GROUPS),
      sample_id = factor(sample_id, targets$sample_id[order(factor(targets$group, FIGURE_GROUPS))])
    ) |>
    ggplot(aes(sample_id, feature, fill = z)) +
    geom_tile() +
    facet_grid(~group, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(
      low = "#2166AC", high = "#B2182B", limits = c(-limit, limit), oob = scales::squish
    ) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_figure() +
    theme(
      axis.text.x = element_text(angle = 90, size = 5, vjust = 0.5),
      axis.text.y = if (row_labels) element_text(size = 5) else element_blank(),
      axis.ticks.y = element_blank(), panel.spacing = unit(1, "mm")
    )
}

# PAS against vehicle per phenotype and sex as a % change with its 95% CI, on a log ratio axis so
# rises and falls of the same size look the same size.
phenotype_effect_plot <- function(effects, title, subtitle) {
  effects |>
    mutate(
      sex = factor(if_else(sex == "F", "Female", "Male"), c("Male", "Female")),
      trait = factor(PHENOTYPE_LABELS[phenotype], rev(PHENOTYPE_LABELS)),
      across(c(pct_change, ci_low, ci_high), \(x) 1 + x / 100)
    ) |>
    ggplot(aes(pct_change, trait, colour = sex)) +
    geom_vline(xintercept = 1, colour = "grey50", linewidth = 0.3) +
    geom_pointrange(aes(xmin = ci_low, xmax = ci_high, shape = P.Value <= ALPHA),
      position = position_dodge(width = 0.6), size = 0.25, linewidth = 0.4, fill = "white"
    ) +
    scale_x_continuous(
      transform = "log2", breaks = c(0.25, 0.5, 1, 2, 4, 8),
      labels = \(x) sprintf("%+.0f%%", 100 * (x - 1))
    ) +
    scale_colour_manual(values = SEX_COLOURS, breaks = c("Female", "Male"), name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c("p > 0.05", "p ≤ 0.05"), name = NULL) +
    labs(x = "PAS vs vehicle, % change (95% CI, log scale)", y = NULL, title = title, subtitle = subtitle) +
    theme_figure()
}

# PCA scores with a shaded 95% t-ellipse per group and a test summary in the top-left corner;
# `arrows` (x, y, label) are drawn from the origin, already scaled to the score axes.
pca_panel <- function(scores, variance, title, subtitle, stats, arrows = NULL) {
  axis <- \(i) sprintf("PC%d (%.1f%%)", i, 100 * variance[i])
  scores <- mutate(scores, group = factor(group, FIGURE_GROUPS))
  plot <- ggplot(scores, aes(PC1, PC2)) +
    stat_ellipse(
      aes(colour = group, fill = group),
      geom = "polygon", type = "t", level = 0.95, linewidth = 0.3, alpha = 0.12
    ) +
    geom_point(aes(colour = group, shape = sex), size = 1.3) +
    annotate(
      "text",
      x = -Inf, y = Inf, label = stats, hjust = -0.05, vjust = 1.3, size = 1.9,
      lineheight = 0.95, colour = "grey20"
    ) +
    scale_colour_manual(values = GROUP_COLOURS, aesthetics = c("colour", "fill"), name = NULL) +
    scale_shape_manual(values = SEX_SHAPES, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.32))) +
    labs(x = axis(1), y = axis(2), title = title, subtitle = subtitle) +
    theme_figure()
  if (is.null(arrows)) {
    return(plot)
  }
  plot +
    geom_segment(
      aes(0, 0, xend = x, yend = y), arrows,
      arrow = grid::arrow(length = grid::unit(1, "mm")), linewidth = 0.3, colour = "grey30"
    ) +
    ggrepel::geom_text_repel(
      aes(x, y, label = label), arrows,
      size = 1.7, colour = "grey20", min.segment.length = Inf, box.padding = 0.15, seed = SEED
    )
}

# One sex's PAS signature (proteins at p <= 0.05, up and down) placed in the other sex's ranking,
# drawn as limma's barcode (Wu & Smyth 2012): members as ticks, the up list above the axis and
# the down list below, with a tricube moving average of membership along the rank. `ranks` holds
# the named moderated t of each target ranking; `signatures` comes from 08_Enrichment.
signature_barcode <- function(signatures, ranks, span = 0.12, floor = 0.34) {
  lane_name <- \(source, ranked_in) sprintf("%s signature\nin %s ranking", str_to_sentence(source), ranked_in)
  lanes <- signatures |>
    mutate(
      lane = factor(lane_name(source, ranked_in), lane_name(c("female", "male"), c("male", "female"))),
      side = if_else(direction == "up", 1, -1),
      ranks = map2(members, ranked_in, \(m, target) {
        order <- rank(-ranks[[target]], ties.method = "first")
        unname(order[str_split_1(m, ";")])
      })
    )
  n_total <- unique(lengths(ranks))
  stopifnot(length(n_total) == 1, all(lengths(lanes$ranks) == lanes$n_ranked))
  ticks <- select(lanes, lane, direction, side, rank = ranks) |>
    unnest(rank)
  worms <- lanes |>
    mutate(worm = map(ranks, \(r) {
      smooth <- limma::tricubeMovingAverage(as.integer(seq_len(n_total) %in% r), span = span)
      tibble(rank = seq_len(n_total), height = floor + (1 - floor) * smooth / max(smooth))
    })) |>
    select(lane, direction, side, worm) |>
    unnest(worm) |>
    mutate(y = side * height, y0 = side * floor)
  labels <- lanes |>
    mutate(
      text = sprintf("%s (n = %d), %.0f%% same way, p %.2g", direction, n_ranked, 100 * same_way, p),
      x = if_else(direction == "up", n_total * 0.995, n_total * 0.005),
      hjust = if_else(direction == "up", 1, 0), y = side * 0.97, vjust = if_else(side > 0, 1, 0)
    )
  ggplot() +
    geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
    geom_segment(
      aes(rank, side * 0.06, xend = rank, yend = side * 0.3, colour = direction), ticks,
      linewidth = 0.2, alpha = 0.8
    ) +
    geom_ribbon(aes(rank, ymin = y0, ymax = y, fill = direction), worms, alpha = 0.25) +
    geom_line(aes(rank, y, colour = direction), worms, linewidth = 0.4) +
    geom_text(
      aes(x, y, label = text, hjust = hjust, vjust = vjust, colour = direction), labels,
      size = 1.9, fontface = "bold", show.legend = FALSE
    ) +
    facet_wrap(~lane, ncol = 1, strip.position = "left") +
    scale_colour_manual(
      values = c(up = DIRECTION_COLOURS[["Up"]], down = DIRECTION_COLOURS[["Down"]]),
      aesthetics = c("colour", "fill"), guide = "none"
    ) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(-1.05, 1.05), breaks = NULL, expand = expansion(0)) +
    theme_figure() +
    theme(
      strip.placement = "outside", strip.text.y.left = element_text(angle = 90, size = 5.5),
      panel.grid.major = element_blank(), panel.spacing.y = unit(2, "mm")
    )
}
