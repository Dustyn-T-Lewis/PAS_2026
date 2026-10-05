# Figure style and the panels more than one figure draws. Sourced by the 06_Figures scripts after
# helpers.R.

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

# p-value histograms, one facet per contrast, with the flat line no effect gives and limma's share
# of true effects (one minus propTrueNull, from 05_Summary/01_Signal).
p_histogram <- function(results, labels, signal, title) {
  hist_data <- results |>
    filter(contrast %in% names(labels)) |>
    mutate(contrast = factor(labels[contrast], labels))
  shares <- signal |>
    filter(term %in% names(labels)) |>
    transmute(
      contrast = factor(labels[term], labels),
      label = sprintf("true effects %.1f%%", 100 * share_true_effects)
    )
  ggplot(hist_data, aes(P.Value)) +
    geom_histogram(breaks = seq(0, 1, 0.05), fill = "grey65", colour = "white", linewidth = 0.2) +
    geom_hline(
      yintercept = nrow(hist_data) / length(labels) * 0.05, linetype = "dashed", linewidth = 0.3
    ) +
    geom_text(
      data = shares, aes(x = 0.97, y = Inf, label = label), hjust = 1, vjust = 1.5, size = 1.9
    ) +
    facet_wrap(~contrast, nrow = 1) +
    scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    labs(x = "p", y = "Proteins", title = title) +
    theme_figure() +
    theme(strip.text = element_text(size = 6, face = "plain"), panel.spacing = unit(3, "mm"))
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

RING_COLLECTIONS <- list(full = c("Hallmark", "Reactome", "KEGG", "GO:BP"), mito = "MitoCarta")
RING_POOLS <- c(
  full = "Hallmark, Reactome (no disease), KEGG reference, GO:BP (no GO Slim)",
  mito = "MitoCarta pathways"
)

study_da <- function(universe) {
  path <- here("02_Differential_Expression", "02_Contrasts", "c_data", paste0("study_", universe))
  suppressMessages(enrichVolcano::read_study(path, species = "Mus musculus"))$da
}

# An enrichVolcano ring of the terms 08_Display chose for one universe and contrast; * marks
# FDR ≤ 0.05, proteins at FDR ≤ 0.10 are coloured in the volcano.
pathway_ring <- function(display, da, universe, contrast, title, subtitle) {
  terms <- filter(display$ring_terms, universe == !!universe, contrast == !!contrast)
  labels <- set_names(
    paste0(
      enrichVolcano::clean_label(terms$term, width = 18), if_else(terms$padj <= ALPHA, " *", "")
    ),
    terms$term
  )
  enrichVolcano::plot_volcano_ring(da, display$display[[universe]],
    contrast = contrast, databases = RING_COLLECTIONS[[universe]], terms = terms$term,
    labels = labels, title = title, subtitle = subtitle, p_threshold = 0.10,
    show_counts = FALSE, score_limits = c(-3, 3), arc_span = "auto", label_size = 1.7,
    label_gap = 1.2, label_headroom = 1.1, axis_size = 1.5,
    theme = enrichVolcano::plot_theme(base_size = 6, base_family = "Helvetica")
  ) +
    theme(
      plot.title = element_text(size = 7.5, face = "bold"),
      plot.subtitle = element_text(size = 5.5, face = "italic", colour = "grey35")
    )
}

# Pathway dot matrix after the HRvLR Figure 1: contrast-coloured column labels, a database square
# beside each term, an optional group strip and an optional term tree, aligned with aplot. `data`
# has term, label, column (factor of contrast labels), NES, padj, database and, if given, group.
pathway_matrix <- function(data, title, subtitle, tree = NULL, group_colours = NULL,
                           group_name = NULL, column_names = NULL) {
  term_order <- if (is.null(tree)) levels(data$term) else tree$labels[tree$order]
  data <- mutate(data, term = factor(term, term_order))
  labels <- set_names(distinct(data, term, label)$label, distinct(data, term, label)$term)
  bands <- distinct(data, column) |> mutate(band = CONTRAST_COLOURS[as.character(column)])
  dots <- ggplot(data, aes(column, term)) +
    geom_tile(
      data = bands, aes(column, y = 1, fill = NULL), fill = bands$band, height = Inf,
      width = 0.95, alpha = 0.15, inherit.aes = FALSE
    ) +
    geom_point(aes(size = -log10(padj), fill = NES), shape = 21, stroke = 0.25) +
    geom_point(
      data = \(d) filter(d, padj <= ALPHA), aes(size = -log10(padj)), shape = 21, fill = NA,
      colour = "black", stroke = 0.6
    ) +
    scale_fill_distiller(
      palette = "PuOr", direction = -1, limits = c(-2.6, 2.6), oob = scales::squish
    ) +
    scale_size_area(max_size = 4, name = "-log10 FDR") +
    scale_y_discrete(position = "right", labels = labels) +
    scale_x_discrete(
      position = "top", expand = expansion(add = 0.6),
      labels = column_names %||% set_names(levels(data$column), levels(data$column))
    ) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_figure() +
    theme(
      axis.text.x.top = element_text(size = 6, face = "bold"),
      axis.ticks = element_blank(),
      axis.text.y = element_text(size = 5.5, colour = "grey15"),
      panel.grid.major.x = element_blank()
    )
  databases <- ggplot(distinct(data, term, database), aes(1, term, fill = database)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = DATABASE_COLOURS, name = "Database") +
    theme_void() +
    theme(legend.key.size = unit(2.5, "mm"), legend.text = element_text(size = 5))
  combined <- aplot::insert_left(dots, databases, width = 0.04)
  if (!is.null(group_colours)) {
    groups <- ggplot(distinct(data, term, group), aes(1, term, fill = group)) +
      geom_tile() +
      scale_fill_manual(values = group_colours, name = group_name) +
      theme_void() +
      theme(legend.key.size = unit(2.5, "mm"), legend.text = element_text(size = 5))
    combined <- aplot::insert_left(combined, groups, width = 0.04)
  }
  if (!is.null(tree)) {
    combined <- aplot::insert_left(
      combined, ggtree::ggtree(ape::as.phylo(tree), linewidth = 0.3),
      width = 0.3
    )
  }
  ggplotify::as.ggplot(combined)
}
