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
DATABASE_COLOURS <- c(
  Hallmark = "#E69F00", Reactome = "#56B4E9", KEGG = "#F0E442", "GO:BP" = "#999999",
  MitoCarta = "#009E73"
)

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
