# Figure style, after HRvLR's: colours sampled from the design panel so a contrast looks the same
# in the schematic and every plot, a minimal theme, and the plots more than one script draws.
# Sourced by helpers.R.

library(ggplot2)
library(patchwork)

SEX_COLOURS <- c(Female = "#D55E00", Male = "#0072B2")
SEX_SHAPES <- c(F = 16, M = 17)
FIGURE_GROUPS <- c("VEH_F", "PAS_F", "VEH_M", "PAS_M")
GROUP_COLOURS <- c(VEH_F = "#FDB863", PAS_F = "#D55E00", VEH_M = "#92C5DE", PAS_M = "#0072B2")
GROUP_LABELS <- c(
  VEH_F = "Female VEH", PAS_F = "Female PAS", VEH_M = "Male VEH", PAS_M = "Male PAS"
)
# Up and down in green and magenta (PiYG), so direction never borrows a group's orange or blue.
DIRECTION_COLOURS <- c(Down = "#C51B7D", Up = "#4D9221")
DIRECTION_LIGHT <- c(Down = "#F1B6DA", Up = "#B8E186")
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
# Each contrast in three tones, as the design panel draws it: line for marks and box borders, fill
# for the box tint, ink for text on white.
CONTRAST_PALETTE <- tibble::tribble(
  ~contrast, ~line, ~fill, ~ink,
  "PAS_vs_VEH_female", "#D55E00", "#FBDCC8", "#A34400",
  "PAS_vs_VEH_male", "#0072B2", "#CCE3F0", "#005282",
  "treatment_by_sex", "#7B3294", "#E8D5EE", "#5A2370",
  "sex_in_VEH", "#737373", "#E5E5E5", "#525252",
  "sex_in_PAS", "#A6A6A6", "#F0F0F0", "#737373"
)
contrast_tone <- \(tone) set_names(CONTRAST_PALETTE[[tone]], CONTRAST_PALETTE$contrast)
DATABASE_COLOURS <- c("GO:BP" = "#999999", Reactome = "#56B4E9", MitoCarta = "#009E73")

# HRvLR's figure theme: minimal, a light major grid, bold title and a plain grey subtitle, with the
# panel letter in a left margin column level with the title. Keys sit under the plot, small, so
# panels in a row share their plot area once aligned.
theme_figure <- function() {
  theme_minimal(base_size = 7, base_family = "Helvetica") +
    theme(
      plot.title = element_text(face = "bold", size = 7.5, margin = margin(b = 1)),
      plot.subtitle = element_text(size = 6, colour = "grey30", margin = margin(b = 3)),
      plot.title.position = "plot",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2),
      legend.position = "bottom",
      legend.key.size = unit(2.5, "mm"),
      legend.text = element_text(size = 6, margin = margin(l = 1, r = 4)),
      legend.margin = margin(0), legend.box.spacing = unit(1.5, "mm"),
      legend.spacing.y = unit(0.5, "mm"),
      plot.tag = element_text(face = "bold", size = 9, vjust = 1, margin = margin(r = 2)),
      plot.tag.location = "margin", plot.tag.position = "left"
    )
}

# Facet strips as the design panel's boxes: the contrast's tint inside a border in its colour, the
# name in bold black; a level outside the palette gets no box.
contrast_strips <- function(levels, side = c("x", "y")) {
  side <- match.arg(side)
  boxes <- ggh4x::elem_list_rect(
    fill = unname(coalesce(contrast_tone("fill")[levels], "transparent")),
    colour = unname(coalesce(contrast_tone("line")[levels], "transparent")),
    linewidth = 0.4
  )
  text <- ggh4x::elem_list_text(colour = "black", face = "bold")
  if (side == "x") {
    ggh4x::strip_themed(background_x = boxes, text_x = text)
  } else {
    ggh4x::strip_themed(background_y = boxes, text_y = text)
  }
}

# A figure panel at the size it takes in its composite, as PDF (cairo) and PNG, into the figure's
# b_reports/panels/; the composite then stitches panels that already look final.
save_panel <- function(plot, figure, name, width_mm, height_mm) {
  dir <- here("04_Figures", figure, "b_reports", "panels")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  save_figure(plot, dir, name, width_mm, height_mm)
}

# Runs one panel script in its own environment and returns what it returns: list(plot, data,
# inputs). Composites call this for each panel, as YvO's stitchers do.
source_panel <- function(figure, name) {
  path <- here("04_Figures", figure, "a_script", "panels", paste0(name, ".R"))
  source(path, local = new.env())$value
}

# PERMANOVA terms as the three lines a PCA panel prints in its corner.
permanova_text <- \(tab) {
  with(filter(tab, term %in% c("sex", "treatment", "sex:treatment")), paste(
    "PERMANOVA",
    paste(sprintf("%s R² %.2f, p = %.2g", sub("sex:treatment", "interaction", term), r2, p),
      collapse = "\n"
    ),
    sep = "\n"
  ))
}

# PCA scores with a shaded 95% t-ellipse per group, a test summary in the top-left corner, and
# arrows (x, y, label) from the origin, already scaled to the score axes.
pca_panel <- function(scores, variance, title, subtitle, stats, arrows) {
  axis <- \(i) sprintf("PC%d (%.1f%%)", i, 100 * variance[i])
  scores <- mutate(scores, group = factor(group, FIGURE_GROUPS))
  ggplot(scores, aes(PC1, PC2)) +
    stat_ellipse(
      aes(colour = group, fill = group),
      geom = "polygon", type = "t", level = 0.95, linewidth = 0.3, alpha = 0.12
    ) +
    geom_point(aes(colour = group, shape = sex), size = 1.3) +
    annotate(
      "label",
      x = -Inf, y = Inf, label = stats, hjust = 0, vjust = 1, size = 1.9, lineheight = 0.95,
      colour = "grey20", fill = alpha("white", 0.75), linewidth = 0,
      label.padding = unit(1, "mm")
    ) +
    scale_colour_manual(
      values = GROUP_COLOURS, labels = GROUP_LABELS, aesthetics = c("colour", "fill"), name = NULL,
      guide = guide_legend(nrow = 2, byrow = TRUE)
    ) +
    scale_shape_manual(values = SEX_SHAPES, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.32))) +
    labs(x = axis(1), y = axis(2), title = title, subtitle = subtitle) +
    theme_figure() +
    geom_segment(
      aes(0, 0, xend = x, yend = y), arrows,
      arrow = grid::arrow(length = grid::unit(1, "mm")), linewidth = 0.3, colour = "grey30"
    ) +
    # The points join the repel layer as empty labels, so the arrow labels move off the mice; a
    # pale backing keeps them legible over the ellipse outlines.
    ggrepel::geom_label_repel(
      aes(x, y, label = label),
      bind_rows(arrows, transmute(scores, x = PC1, y = PC2, label = "")),
      size = 1.7, colour = "grey20", fill = alpha("white", 0.75), linewidth = 0,
      label.padding = 0.08, min.segment.length = Inf, box.padding = 0.2, point.padding = 0.1,
      max.overlaps = Inf, seed = SEED
    )
}

# One sex's PAS signature (proteins at p <= 0.05, up and down) placed in the other sex's ranking,
# drawn as limma's barcode (Wu & Smyth 2012): members as ticks, the up list above the axis and
# the down list below, with a tricube moving average of membership along the rank. ranks holds
# the named moderated t of each target ranking; signatures comes from 02_Enrichment.
signature_barcode <- function(signatures, ranks, span = 0.12, min_height = 0.34) {
  lane_name <- \(source, ranked_in) {
    sprintf("%s signature\nin %s ranking", str_to_sentence(source), ranked_in)
  }
  lanes <- signatures |>
    mutate(
      lane = factor(
        lane_name(source, ranked_in), lane_name(c("female", "male"), c("male", "female"))
      ),
      side = if_else(direction == "up", 1, -1),
      ranks = map2(members, ranked_in, \(m, target) {
        position <- rank(-ranks[[target]], ties.method = "first")
        unname(position[str_split_1(m, ";")])
      })
    )
  n_total <- unique(lengths(ranks))
  stopifnot(length(n_total) == 1, all(lengths(lanes$ranks) == lanes$n_ranked))
  ticks <- select(lanes, lane, direction, side, rank = ranks) |>
    unnest(rank)
  worms <- lanes |>
    mutate(worm = map(ranks, \(r) {
      smooth <- limma::tricubeMovingAverage(as.integer(seq_len(n_total) %in% r), span = span)
      tibble(rank = seq_len(n_total), height = min_height + (1 - min_height) * smooth / max(smooth))
    })) |>
    select(lane, direction, side, worm) |>
    unnest(worm) |>
    mutate(y = side * height, y0 = side * min_height)
  labels <- lanes |>
    mutate(
      text = sprintf(
        "%s (n = %d), %.0f%% same way, fry p %.2g", direction, n_ranked, 100 * same_way, p
      ),
      x = if_else(direction == "up", n_total * 0.995, n_total * 0.005),
      hjust = if_else(direction == "up", 1, 0), y = side * 1.38, vjust = if_else(side > 0, 1, 0)
    )
  ggplot() +
    geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
    geom_segment(
      aes(rank, side * 0.06, xend = rank, yend = side * 0.3, colour = direction), ticks,
      linewidth = 0.2, alpha = 0.8, show.legend = FALSE
    ) +
    geom_ribbon(
      aes(rank, ymin = y0, ymax = y, fill = direction, colour = direction), worms,
      alpha = 0.25, linewidth = 0, key_glyph = "rect"
    ) +
    geom_line(aes(rank, y, colour = direction), worms, linewidth = 0.4, show.legend = FALSE) +
    geom_text(
      aes(x, y, label = text, hjust = hjust, vjust = vjust, colour = direction), labels,
      size = 1.9, fontface = "bold", show.legend = FALSE
    ) +
    facet_wrap(~lane, ncol = 1, strip.position = "left") +
    scale_colour_manual(
      values = c(up = DIRECTION_COLOURS[["Up"]], down = DIRECTION_COLOURS[["Down"]]),
      labels = c(up = "Up in its own sex", down = "Down in its own sex"),
      aesthetics = c("colour", "fill"), name = NULL,
      guide = guide_legend(nrow = 2, reverse = TRUE)
    ) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(-1.42, 1.42), breaks = NULL, expand = expansion(0)) +
    theme_figure() +
    theme(
      # Inside, not outside: patchwork then puts this axis in the same row as its neighbours'.
      strip.placement = "inside", strip.text.y.left = element_text(angle = 90, size = 5.5),
      panel.grid.major = element_blank(), panel.spacing.y = unit(2, "mm")
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

# Boxed contrast labels, as the design panel draws them: the contrast's tint inside a border in its
# colour, the text in bold black, the box only as large as its text. data holds x, y, label and
# contrast; extra arguments (hjust, vjust) pass to geom_label.
contrast_labels <- function(data, ...) {
  geom_label(
    aes(x, y,
      label = label, fill = I(contrast_tone("fill")[contrast]),
      colour = I(contrast_tone("line")[contrast])
    ),
    data,
    text.colour = "black", fontface = "bold", size = 1.8, linewidth = 0.4,
    label.r = unit(0.5, "mm"), label.padding = unit(0.45, "mm"), lineheight = 0.85,
    inherit.aes = FALSE, show.legend = FALSE, ...
  )
}

# A step's preview PDF: one landscape page per named function, each drawing with a package's own
# plot or base graphics; cairo, so Π and ≤ render. Returns the rows of the workbook's plot_index.
preview_pdf <- function(file, pages) {
  cairo_pdf(file, width = 7, height = 5, onefile = TRUE)
  on.exit(dev.off())
  defaults <- par(no.readonly = TRUE)
  walk(pages, \(draw) {
    draw()
    par(defaults)
  })
  tibble(file = basename(file), page = seq_along(pages), shows = names(pages))
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
