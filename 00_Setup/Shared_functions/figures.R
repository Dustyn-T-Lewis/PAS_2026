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
GROUP_SHORT <- c(VEH_F = "F VEH", PAS_F = "F PAS", VEH_M = "M VEH", PAS_M = "M PAS")
# Up and down in green and magenta (PiYG), so direction never borrows a group's orange or blue.
DIRECTION_COLOURS <- c(Down = "#C51B7D", Up = "#4D9221")
DIRECTION_LIGHT <- c(Down = "#F1B6DA", Up = "#B8E186")
PHENOTYPE_LABELS <- c(
  csa_mean = "CSA", capillary_density = "Capillary density",
  satellite_cells_per_fiber = "Satellite cells / fibre",
  conductance_sol_pyr = "Soleus pyruvate", conductance_pla_pyr = "Plantaris pyruvate",
  conductance_sol_oct = "Soleus octanoyl"
)
TRAIT_SHORT <- c(
  csa_mean = "CSA", capillary_density = "Capillaries", satellite_cells_per_fiber = "Sat. cells",
  conductance_sol_pyr = "Sol. pyr", conductance_pla_pyr = "Pla. pyr",
  conductance_sol_oct = "Sol. oct"
)
MODEL_LABELS <- c(pooled = "Pooled", treatment = "Within treatment")
CONTRAST_LABELS <- c(
  PAS_vs_VEH_female = "PAS vs VEH, female", PAS_vs_VEH_male = "PAS vs VEH, male",
  treatment_by_sex = "Treatment × sex", sex_in_VEH = "Sex in VEH", sex_in_PAS = "Sex in PAS"
)
CONTRAST_SHORT <- c(
  PAS_vs_VEH_female = "Female", PAS_vs_VEH_male = "Male", treatment_by_sex = "F − M"
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
DATABASE_COLOURS <- c(
  "GO:BP" = "#999999", Reactome = "#56B4E9", Hallmark = "#E69F00", KEGG = "#CC79A7",
  MitoCarta = "#009E73"
)

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

# Dot size and ring shared by the dot matrices: area grows with -log10 FDR from a visible minimum,
# and a thick ring marks fgsea FDR <= 0.05.
fdr_scales <- function(fdr_cap) {
  list(
    scale_size(
      range = c(0.5, 3.2), limits = c(0, fdr_cap), breaks = c(1, 2, fdr_cap),
      labels = c("1", "2", paste0(fdr_cap, "+")), name = "-log10 FDR"
    ),
    scale_discrete_manual(
      "stroke",
      values = c(`TRUE` = 0.7, `FALSE` = 0.15), breaks = "TRUE", labels = "FDR ≤ 0.05",
      name = NULL, guide = guide_legend(override.aes = list(size = 2.4, fill = "white"))
    )
  )
}

# Pathway sets as a dot matrix, after F01's old panel f: NES as fill (magenta down, green up),
# -log10 FDR as size, the collection as a square before the first column. cells holds term,
# column, col_group, NES and fdr; rows holds term, label, collection and row_group, one row per
# set and group in display order. strip styles the row and column strips (ggh4x::strip_themed).
# Every collection is listed in the key and nes_limit can be fixed, so two matrices can share
# one key.
pathway_matrix <- function(cells, rows, title, subtitle, strip = ggh4x::strip_vanilla(),
                           fdr_cap = 4, nes_limit = NULL) {
  rows <- mutate(rows, row = factor(paste(row_group, term), rev(paste(row_group, term))))
  cells <- cells |>
    inner_join(select(rows, term, row, row_group), by = "term", relationship = "many-to-many") |>
    mutate(log_fdr = pmin(-log10(fdr), fdr_cap))
  squares <- mutate(rows, col_group = factor(levels(cells$col_group)[1], levels(cells$col_group)))
  nes_limit <- nes_limit %||% (ceiling(max(abs(cells$NES)) * 10) / 10)
  ggplot(cells, aes(column, row)) +
    geom_point(
      aes(size = log_fdr, fill = NES, stroke = fdr <= ALPHA),
      shape = 21, colour = "grey15"
    ) +
    geom_point(aes(x = 0.3, colour = collection), squares, shape = 15, size = 1.5) +
    # Every collection once, invisible, so each has a key even where a matrix lacks it and two
    # matrices' keys match.
    geom_point(
      aes(x = 0.3, y = row, colour = collection),
      tibble(
        row = rows$row[1], row_group = rows$row_group[1], col_group = squares$col_group[1],
        collection = names(DATABASE_COLOURS)
      ),
      alpha = 0, inherit.aes = FALSE
    ) +
    ggh4x::facet_grid2(row_group ~ col_group,
      scales = "free", space = "free", switch = "y", strip = strip
    ) +
    scale_y_discrete(labels = deframe(select(rows, row, label))) +
    scale_x_discrete(expand = expansion(add = c(0.8, 0.5))) +
    scale_colour_manual(
      values = DATABASE_COLOURS, limits = names(DATABASE_COLOURS), drop = FALSE, name = NULL
    ) +
    scale_fill_gradient2(
      low = DIRECTION_COLOURS[["Down"]], high = DIRECTION_COLOURS[["Up"]],
      limits = c(-nes_limit, nes_limit), breaks = c(-2, 0, 2), name = "NES"
    ) +
    fdr_scales(fdr_cap) +
    guides(
      fill = guide_colourbar(order = 1, barwidth = unit(14, "mm"), barheight = unit(2, "mm")),
      size = guide_legend(order = 2),
      colour = guide_legend(
        order = 3, nrow = 2,
        override.aes = list(shape = 15, size = 1.5, alpha = 1)
      )
    ) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_figure() +
    theme(
      axis.text.y = element_text(size = 5.2), axis.text.x = element_text(size = 5.5),
      strip.placement = "outside", strip.text.y.left = element_text(angle = 90),
      strip.text = element_text(face = "bold", size = 5.5), panel.spacing = unit(0.8, "mm"),
      legend.box = "vertical", legend.box.just = "left"
    )
}

# Display names for pathway sets: the collection prefix dropped, sentence case, shortened.
set_label <- \(term, width = 48) str_trunc(enrichVolcano::clean_label(term, width = 200), width)

# Proteins by mice as a ComplexHeatmap, columns split by group under a group-colour bar, values
# already centred (log2 against a reference), magenta below and green above. Returns a grob
# patchwork can place with wrap_elements(); right is an optional rowAnnotation.
protein_heatmap <- function(values, groups, row_split = NULL, right = NULL, limit = 1,
                            name = "log2 vs\nown-sex VEH", width_mm = 120, height_mm = 120,
                            legend = TRUE) {
  groups <- factor(groups, FIGURE_GROUPS)
  colours <- circlize::colorRamp2(
    c(-limit, 0, limit), c(DIRECTION_COLOURS[["Down"]], "white", DIRECTION_COLOURS[["Up"]])
  )
  top <- ComplexHeatmap::HeatmapAnnotation(
    group = groups, col = list(group = GROUP_COLOURS), show_annotation_name = FALSE,
    show_legend = FALSE, simple_anno_size = unit(1.6, "mm")
  )
  small <- grid::gpar(fontsize = 5.5)
  heatmap <- ComplexHeatmap::Heatmap(values,
    name = name, col = colours, top_annotation = top, right_annotation = right,
    column_split = groups, cluster_column_slices = FALSE, column_title_gp = small,
    column_title = GROUP_SHORT[levels(droplevels(groups))],
    column_labels = rep("", ncol(values)), row_split = row_split, cluster_row_slices = FALSE,
    row_title_gp = small, row_title_rot = 0, show_row_names = FALSE, show_row_dend = FALSE,
    show_column_dend = FALSE, column_gap = unit(0.6, "mm"), row_gap = unit(0.8, "mm"),
    show_heatmap_legend = legend,
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = 5.5, fontface = "bold"), labels_gp = small,
      grid_width = unit(2, "mm")
    )
  )
  grid::grid.grabExpr(
    ComplexHeatmap::draw(heatmap,
      heatmap_legend_side = "right", annotation_legend_side = "right",
      padding = unit(c(1, 1, 4, 1), "mm")
    ),
    width = width_mm / 25.4, height = height_mm / 25.4
  )
}

# Each comparison contrast's sets at FDR <= 0.05 that collapsePathways keeps, listed once under the
# first contrast that selects them and drawn in all three columns. results is 02_Enrichment's
# table, already filtered to the collections shown.
contrast_pathway_matrix <- function(results, title, subtitle, nes_limit = NULL) {
  rows <- results |>
    filter(padj <= ALPHA, dedup_status == "kept") |>
    arrange(factor(contrast, COMPARISON_CONTRASTS), p) |>
    distinct(term, .keep_all = TRUE) |>
    transmute(
      term, collection,
      label = set_label(term),
      row_group = factor(CONTRAST_SHORT[contrast], CONTRAST_SHORT)
    )
  rows$row_group <- droplevels(rows$row_group)
  cells <- results |>
    filter(term %in% rows$term) |>
    transmute(
      term, NES,
      fdr = padj, column = factor(CONTRAST_SHORT[contrast], CONTRAST_SHORT),
      col_group = factor("")
    )
  present <- names(CONTRAST_SHORT)[CONTRAST_SHORT %in% levels(rows$row_group)]
  plot <- pathway_matrix(cells, rows, title, subtitle,
    strip = contrast_strips(present, "y"), nes_limit = nes_limit
  ) +
    theme(strip.text.x = element_blank())
  list(plot = plot, rows = rows, cells = cells)
}

# One sex's pathway associations: phenotype by phenotype, the n sets with the smallest fgsea FDR,
# among those at FDR <= 0.05 that collapsePathways keeps under either model, that no earlier
# phenotype already listed; each set drawn across every phenotype under both models. assoc is
# 02_Enrichment's association_fgsea.
phenotype_pathway_matrix <- function(assoc, sex, title, subtitle, n = 4, nes_limit = NULL) {
  assoc <- filter(assoc, sex == !!sex, model %in% names(MODEL_LABELS))
  candidates <- assoc |>
    filter(fgsea_fdr <= ALPHA, dedup_status == "kept") |>
    summarise(best = min(fgsea_fdr), best_p = min(fgsea_p), .by = c(phenotype, term, collection)) |>
    arrange(best, best_p)
  rows <- reduce(names(TRAIT_SHORT), \(listed, trait) {
    candidates |>
      filter(phenotype == trait, !term %in% listed$term) |>
      slice_head(n = n) |>
      bind_rows(x = listed)
  }, .init = candidates[0, ]) |>
    transmute(
      term, collection,
      label = set_label(term, 44),
      row_group = factor(TRAIT_SHORT[phenotype], TRAIT_SHORT)
    )
  cells <- assoc |>
    filter(term %in% rows$term) |>
    transmute(
      term, NES,
      fdr = fgsea_fdr, column = factor(TRAIT_SHORT[phenotype], TRAIT_SHORT),
      col_group = factor(MODEL_LABELS[model], MODEL_LABELS)
    )
  plot <- pathway_matrix(cells, rows, title, subtitle, nes_limit = nes_limit) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  list(plot = plot, rows = rows, cells = cells)
}

# A figure's panel scripts, run in order, each plot carrying its letter.
run_panels <- function(figure, panels) {
  map(set_names(panels), \(name) {
    panel <- source_panel(figure, name)
    panel$plot <- panel$plot + labs(tag = letters[match(name, panels)])
    panel
  })
}

# A composite's workbook: one sheet per panel, then every input the panels read.
write_figure_workbook <- function(panels, holds, file) {
  inputs <- unlist(unname(map(panels, "inputs")))
  write_workbook(map(panels, "data"), holds, file, inputs[!duplicated(inputs)])
}

# One contrast's proteins at Π <= 0.05 (full proteome) across all twenty mice, each centred on its
# own sex's vehicle mean, split by direction (for the interaction, by which sex's PAS effect is
# larger). fit is 01_Fit's output.
pi_heatmap_panel <- function(fit, contrast, title, width_mm = 58, height_mm = 120,
                             legend = TRUE) {
  selected <- fit$protein_contrasts |>
    filter(universe == "full", contrast == !!contrast, pi_score <= ALPHA) |>
    mutate(direction = if_else(logFC > 0, 1, 2)) |>
    arrange(direction)
  targets <- arrange(fit$targets, factor(group, FIGURE_GROUPS))
  values <- fit$E_norm$full[selected$feature, targets$sample_id]
  vehicle_mean <- sapply(c(F = "VEH_F", M = "VEH_M"), \(g) rowMeans(values[, targets$group == g]))
  names <- if (contrast == "treatment_by_sex") c("F > M", "F < M") else c("up", "down")
  selected <- mutate(selected, direction = factor(names[direction], names))
  counts <- count(selected, direction)
  split <- factor(
    sprintf("%s (%d)", selected$direction, counts$n[match(selected$direction, counts$direction)]),
    sprintf("%s (%d)", counts$direction, counts$n)
  )
  heatmap <- protein_heatmap(values - vehicle_mean[, targets$sex], targets$group,
    row_split = split, limit = 1.5, width_mm = width_mm, height_mm = height_mm, legend = legend
  )
  plot <- wrap_elements(plot = heatmap) +
    labs(title = title, subtitle = sprintf("Π ≤ 0.05: %d proteins", nrow(selected))) +
    theme_figure()
  list(
    plot = plot,
    data = select(selected, feature, gene, direction, logFC, P.Value, adj.P.Val, pi_score)
  )
}

# Sets by contrast as dots without a collection column, for side lists: NES as fill, -log10 FDR
# as size, rows grouped by row_group. cells holds row (display order), label, row_group, contrast,
# NES and fdr.
set_dot_matrix <- function(cells, title, subtitle, fdr_cap = 4, nes_limit = 2.5) {
  cells |>
    mutate(
      column = factor(CONTRAST_SHORT[contrast], CONTRAST_SHORT),
      log_fdr = pmin(-log10(fdr), fdr_cap)
    ) |>
    ggplot(aes(column, row)) +
    geom_point(
      aes(size = log_fdr, fill = NES, stroke = fdr <= ALPHA),
      shape = 21, colour = "grey15"
    ) +
    facet_grid(row_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_discrete(labels = deframe(distinct(cells, row, label))) +
    scale_fill_gradient2(
      low = DIRECTION_COLOURS[["Down"]], high = DIRECTION_COLOURS[["Up"]],
      limits = c(-nes_limit, nes_limit), breaks = c(-2, 0, 2), name = "NES",
      guide = guide_colourbar(barwidth = unit(14, "mm"), barheight = unit(2, "mm"))
    ) +
    fdr_scales(fdr_cap) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_figure() +
    theme(
      axis.text.y = element_text(size = 5.2), strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 5.2),
      panel.spacing = unit(0.8, "mm"), legend.box = "vertical"
    )
}
