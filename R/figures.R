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

# pipeline.png: one box per stage folder listing its steps (00_Input lists the inputs read), and
# an arrow wherever one stage reads another's saved files, taken from the input_manifest sheet of
# every step's workbook. A stage sits one column right of the furthest stage it reads; an arrow
# that skips a column runs below the boxes it passes.
draw_pipeline <- function(steps, file) {
  step_dirs <- dirname(dirname(steps))
  stage_of <- \(path) str_extract(path, "^[^/]+")
  reads <- map(step_dirs, \(d) {
    workbook <- here(d, "c_data", paste0(tolower(basename(d)), ".xlsx"))
    tibble(stage = stage_of(d), path = readxl::read_excel(workbook, "input_manifest")$path)
  }) |>
    list_rbind()
  edges <- reads |>
    distinct(from = stage_of(path), to = stage) |>
    filter(from != to)
  inputs <- sort(unique(na.omit(str_extract(reads$path, "(?<=^00_Input/)[^/]+"))))
  items <- tibble(stage = c(rep("00_Input", length(inputs)), stage_of(step_dirs))) |>
    mutate(item = c(inputs, basename(step_dirs)))

  # Steps run in dependency order, so every stage's sources already have a column.
  stages <- unique(items$stage)
  column <- set_names(rep(0, length(stages)), stages)
  for (s in stages) {
    sources <- edges$from[edges$to == s]
    if (length(sources) > 0) column[[s]] <- 1 + max(column[sources])
  }

  width <- 2.6
  step_x <- 4
  gap_y <- 0.5
  boxes <- items |>
    summarise(text = paste(item, collapse = "\n"), n = n(), .by = stage) |>
    mutate(column = column[stage], x = column * step_x, height = 0.5 + 0.22 * n) |>
    mutate(
      top = -lag(cumsum(height + gap_y), default = 0),
      bottom = top - height,
      shift = (max(top) + min(bottom)) / 2,
      top = top - shift,
      bottom = bottom - shift,
      .by = column
    )

  lane_below <- \(from_col, to_col) min(boxes$bottom[between(boxes$column, from_col, to_col - 1)])
  paths <- edges |>
    left_join(select(boxes, from = stage, fx = x, fbottom = bottom, ftop = top, fcol = column),
      by = "from"
    ) |>
    left_join(select(boxes, to = stage, tx = x, ttop = top, tbottom = bottom, tcol = column),
      by = "to"
    ) |>
    mutate(
      fy = (ftop + fbottom) / 2, ty = (ttop + tbottom) / 2, skip = tcol - fcol > 1,
      bus = tx - width / 2 - 0.3 * dense_rank(desc(fcol)), .by = tcol
    ) |>
    mutate(lane = map2_dbl(fcol, tcol, lane_below) - 0.3 * dense_rank(fcol), .by = skip) |>
    mutate(id = row_number()) |>
    pmap(\(id, from, fx, fy, fbottom, tx, ty, bus, lane, skip, ...) {
      start <- if (skip) list(c(fx, fx), c(fbottom, lane)) else list(fx + width / 2, fy)
      tibble(
        id = id, from = from,
        x = c(start[[1]], bus, bus, tx - width / 2),
        y = c(start[[2]], if (skip) lane else fy, ty, ty)
      )
    }) |>
    list_rbind()

  colours <- set_names(grDevices::palette.colors(length(stages) + 1, "Okabe-Ito")[-1], stages)
  figure <- ggplot() +
    geom_path(aes(x, y, group = id, colour = from),
      data = paths, linewidth = 0.45,
      arrow = arrow(length = unit(1.6, "mm"), type = "closed")
    ) +
    geom_rect(aes(
      xmin = x - width / 2, xmax = x + width / 2, ymin = bottom, ymax = top,
      colour = stage
    ), data = boxes, fill = "white", linewidth = 0.7) +
    geom_text(aes(x - width / 2 + 0.12, top - 0.12, label = stage),
      data = boxes, hjust = 0, vjust = 1, size = 3, fontface = "bold"
    ) +
    geom_text(aes(x - width / 2 + 0.12, top - 0.45, label = text),
      data = boxes, hjust = 0, vjust = 1, size = 2.4, lineheight = 0.95, colour = "grey25"
    ) +
    scale_colour_manual(values = colours, guide = "none") +
    labs(
      title = "Pipeline: Rscript run_all.R",
      caption = paste(
        "Arrows: a stage reads saved files of another (input_manifest sheet of each workbook).",
        "Each step folder holds a_script/ (script), b_reports/ (PDF) and c_data/ (workbook, .rds).",
        sep = "\n"
      )
    ) +
    coord_equal(clip = "off") +
    theme_void(base_size = 9, base_family = "Helvetica") +
    theme(
      plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(8, 8, 8, 8)
    )
  ggsave(file, figure, width = 250, height = 140, units = "mm", dpi = 300)
}
