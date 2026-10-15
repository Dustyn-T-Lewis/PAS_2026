# Figure style and the plots more than one script draws, for figures and step reports alike.
# Sourced by R/helpers.R.

library(ggplot2)
library(patchwork)

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

# Framed panels with a light grid, bold title and italic subtitle.
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

# A finished figure as PDF (cairo, which carries Π and ≤) and as a 600 dpi PNG.
save_figure <- function(figure, dir, name, width_mm, height_mm) {
  save <- \(ext, ...) {
    ggsave(file.path(dir, paste0(name, ext)), figure,
      width = width_mm, height = height_mm, units = "mm", ...
    )
  }
  save(".pdf", device = grDevices::cairo_pdf)
  save(".png", dpi = 600)
}

# A step report: one landscape page per named ggplot, cairo so Π and ≤ render. Returns the
# rows of the workbook's plot_index sheet.
write_report <- function(file, pages) {
  cairo_pdf(file, width = 7, height = 5, onefile = TRUE)
  on.exit(dev.off())
  walk(pages, print)
  tibble(file = basename(file), page = seq_along(pages), shows = names(pages))
}

# A tint of a colour, moved the given share of the way towards white.
lighten <- function(colour, amount) {
  grDevices::colorRampPalette(c(colour, "white"))(101)[round(100 * amount) + 1]
}

# Rows ordered by hierarchical clustering on 1 - Pearson r, average linkage: display only.
cluster_order <- \(m) hclust(as.dist(1 - cor(t(m), use = "pairwise.complete.obs")), "average")$order

# Rows z-scored, ordered by cluster_order(), mice in group order.
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
      low = DIRECTION_COLOURS[["Down"]], high = DIRECTION_COLOURS[["Up"]],
      limits = c(-limit, limit), oob = scales::squish
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
    scale_shape_manual(
      values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "p ≤ 0.05", `FALSE` = "p > 0.05"),
      name = NULL
    ) +
    labs(
      x = "PAS vs vehicle, % change (95% CI, log scale)", y = NULL, title = title,
      subtitle = subtitle
    ) +
    theme_figure()
}
