# S3's columns, one per comparison contrast (a volcano over the per-mouse heatmap and its cluster
# ORA bars), and the key they share; sourced by the S3 panel scripts and S3.R.

# The fill reaches full magenta and green at this log2 distance from the own-sex vehicle mean, in
# the heatmaps and in the key.
HEATMAP_LIMIT <- 1
heatmap_colours <- circlize::colorRamp2(
  c(-HEATMAP_LIMIT, 0, HEATMAP_LIMIT),
  c(DIRECTION_COLOURS[["Down"]], "white", DIRECTION_COLOURS[["Up"]])
)
VOLCANO_STATUS <- c(
  down = paste("Down, Π ≤", ALPHA), up = paste("Up, Π ≤", ALPHA), none = paste("Π >", ALPHA)
)
VOLCANO_COLOURS <- set_names(
  c(DIRECTION_COLOURS[["Down"]], DIRECTION_COLOURS[["Up"]], "grey80"), VOLCANO_STATUS
)
HEATMAP_MM <- 112

# One contrast's column. On top, every protein in the full proteome as log2 FC against -log10 p,
# Π <= 0.05 coloured by direction and FDR <= 0.05 ringed and named, on the axes all three
# comparison contrasts share. Below, its proteins at Π <= 0.05 across every mouse, each centred
# on its own sex's vehicle mean, grey where not detected, from 02_Enrichment's heatmap_clusters:
# one row slice per co-expression cluster in drawing order, the proteins the cut left unassigned
# under their direction's last cluster, and the direction (for the interaction, which sex's PAS
# effect is larger) as a thin strip on the left. Columns split by group and clustered within it.
# Beside each cluster, in anno_zoom()'s box, the n_sets sets with the smallest ORA p (then FDR)
# among all it was tested against (cluster_ora), whether or not they reach p <= ALPHA, as -log10 p
# bars on the cluster's own range, filled by database and outlined at FDR <= ALPHA. fit is
# 01_Fit's output, quant 02_Quantification's and enrichment 02_Enrichment's. width_mm and
# height_mm are the column's size, returned for save_figure; the bars take zoom_mm of the width.
pi_heatmap_panel <- function(fit, quant, enrichment, contrast, width_mm = 80, height_mm = 172,
                             zoom_mm = 41, n_sets = 3) {
  proteins <- fit$protein_contrasts |>
    filter(universe == "full", contrast %in% COMPARISON_CONTRASTS) |>
    mutate(
      status = case_when(
        pi_score > ALPHA ~ VOLCANO_STATUS[["none"]], logFC > 0 ~ VOLCANO_STATUS[["up"]],
        .default = VOLCANO_STATUS[["down"]]
      ),
      status = factor(status, VOLCANO_STATUS)
    )
  limits <- list(x = range(proteins$logFC), y = c(0, max(-log10(proteins$P.Value))))
  proteins <- proteins |>
    filter(contrast == !!contrast) |>
    arrange(desc(status))
  named <- filter(proteins, adj.P.Val <= ALPHA)

  names <- if (contrast == "treatment_by_sex") c("F > M", "F < M") else c("up", "down")
  targets <- arrange(fit$targets, factor(group, FIGURE_GROUPS))
  selected <- enrichment$heatmap_clusters |>
    filter(contrast == !!contrast) |>
    arrange(row) |>
    mutate(direction = factor(names[match(direction, c("up", "down"))], names), slice = cluster) |>
    fill(slice, .by = direction) |>
    arrange(factor(slice, unique(slice)))
  slices <- factor(selected$slice, unique(selected$slice))
  values <- fit$E_norm$full[selected$feature, targets$sample_id]
  centred <- centre_on_vehicle(values, targets)
  detected <- quant$proteins$other$n.observations[selected$feature, targets$sample_id] > 0
  top <- enrichment$cluster_ora |>
    filter(contrast == !!contrast) |>
    arrange(p, fdr) |>
    slice_head(n = n_sets, by = cluster) |>
    mutate(label = sprintf("%s (%d/%d)", set_label(term, 36), overlap, size))

  volcano <- ggplot(proteins, aes(logFC, -log10(P.Value))) +
    geom_point(aes(colour = status), size = 0.5, alpha = 0.8) +
    geom_point(data = named, shape = 21, size = 1.6, stroke = 0.4) +
    ggrepel::geom_text_repel(aes(label = gene), named,
      size = 5 / .pt, min.segment.length = 0, segment.size = 0.2, point.padding = 0.4, seed = SEED
    ) +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", linewidth = 0.25) +
    scale_colour_manual(values = VOLCANO_COLOURS, guide = "none") +
    coord_cartesian(xlim = limits$x, ylim = limits$y) +
    labs(
      x = "log2 fold change", y = "-log10 p", title = CONTRAST_LABELS[[contrast]],
      subtitle = sprintf(
        "Π ≤ %s: %d proteins; dashed: p = %s; ringed: FDR ≤ %s", ALPHA, nrow(selected), ALPHA, ALPHA
      )
    ) +
    theme_figure()

  small <- grid::gpar(fontsize = 5.5)
  line_mm <- 1.8
  bar_mm <- 1.5
  slot_mm <- line_mm + bar_mm + 0.7
  axis_mm <- 4
  box_mm <- as.vector(table(factor(top$cluster, levels(slices)))) * slot_mm + axis_mm
  # Each cluster's sets from the top of its box, the name above its bar, on the cluster's own
  # -log10 p range (always past p = ALPHA, dotted on each bar), with a small axis at the foot.
  bars <- function(index, nm) {
    sets <- filter(top, cluster == nm)
    breaks <- pretty(c(0, -log10(c(sets$p, ALPHA))), 2)
    reach <- max(breaks)
    grid::pushViewport(grid::viewport(
      x = unit(0.8, "mm"), width = unit(1, "npc") - unit(2.3, "mm"), just = "left",
      xscale = c(0, reach)
    ))
    text_y <- unit(1, "npc") - unit((seq_len(nrow(sets)) - 1) * slot_mm, "mm")
    bar_y <- text_y - unit(line_mm, "mm")
    grid::grid.text(sets$label,
      x = 0, y = text_y, just = c("left", "top"), gp = grid::gpar(fontsize = 5)
    )
    grid::grid.rect(
      x = unit(0, "native"), y = bar_y, width = unit(-log10(sets$p), "native"),
      height = unit(bar_mm, "mm"), just = c("left", "top"),
      gp = grid::gpar(
        fill = DATABASE_COLOURS[sets$collection], lwd = 0.7,
        col = if_else(sets$fdr <= ALPHA, "black", NA_character_)
      )
    )
    grid::grid.segments(
      unit(-log10(ALPHA), "native"), bar_y + unit(0.3, "mm"),
      unit(-log10(ALPHA), "native"), bar_y - unit(bar_mm + 0.3, "mm"),
      gp = grid::gpar(lty = "11", lwd = 0.7)
    )
    grid::pushViewport(grid::viewport(
      y = unit(axis_mm - 0.6, "mm"), height = unit(0, "mm"), xscale = c(0, reach)
    ))
    grid::grid.xaxis(
      at = breaks,
      gp = grid::gpar(fontsize = 4.5, lwd = 0.4, col = "grey30"),
      edits = grid::gEdit("ticks", y1 = unit(-0.5, "mm"))
    )
    grid::popViewport(2)
  }

  block_of <- selected$direction[!duplicated(slices)]
  cluster_n <- table(factor(selected$cluster, levels(slices)))
  groups <- factor(targets$group, FIGURE_GROUPS)
  not_detected <- \(j, i, x, y, w, h, fill) {
    hide <- !ComplexHeatmap::pindex(detected, i, j)
    if (any(hide)) {
      grid::grid.rect(x[hide], y[hide], w[hide], h[hide],
        gp = grid::gpar(fill = "grey75", col = NA)
      )
    }
  }
  heatmap <- ComplexHeatmap::Heatmap(centred,
    col = heatmap_colours, column_split = groups, cluster_column_slices = FALSE,
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      group = groups, col = list(group = GROUP_COLOURS), show_annotation_name = FALSE,
      show_legend = FALSE, simple_anno_size = unit(1.6, "mm")
    ),
    left_annotation = ComplexHeatmap::rowAnnotation(
      direction = selected$direction,
      col = list(direction = set_names(DIRECTION_COLOURS[c("Up", "Down")], names)),
      show_annotation_name = FALSE, show_legend = FALSE, simple_anno_size = unit(1.2, "mm")
    ),
    right_annotation = ComplexHeatmap::rowAnnotation(
      ora = ComplexHeatmap::anno_zoom(
        align_to = split(seq_along(slices), slices), which = "row", panel_fun = bars,
        size = unit(box_mm, "mm"), gap = unit(1.2, "mm"), width = unit(zoom_mm, "mm"),
        link_width = unit(2.5, "mm"), link_gp = grid::gpar(col = "grey55", fill = NA, lwd = 0.5)
      ),
      show_annotation_name = FALSE
    ),
    column_title = GROUP_LABELS[levels(droplevels(groups))], column_title_gp = small,
    column_title_rot = 90, show_column_names = FALSE, show_column_dend = FALSE,
    row_split = slices, cluster_rows = FALSE, show_row_names = FALSE,
    row_title = sprintf("%s (%d)", levels(slices), cluster_n), row_title_rot = 0,
    row_title_gp = grid::gpar(fontsize = 5),
    row_gap = unit(if_else(head(block_of, -1) == tail(block_of, -1), 1, 2.5), "mm"),
    column_gap = unit(0.6, "mm"), layer_fun = not_detected, show_heatmap_legend = FALSE
  )
  # grid.grabExpr returns to the device open before it, and with none open R starts Rplots.pdf.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  grob <- grid::grid.grabExpr(
    ComplexHeatmap::draw(heatmap, padding = unit(c(1, 1, 1, 1), "mm")),
    width = width_mm / 25.4, height = HEATMAP_MM / 25.4
  )
  # The heatmap keeps the size it was drawn at; the volcano takes the rest of the column.
  column <- volcano / wrap_elements(full = grob) +
    plot_layout(heights = unit(c(1, HEATMAP_MM), c("null", "mm")))
  cells <- as_tibble(centred, rownames = "feature") |>
    pivot_longer(-feature, names_to = "sample_id", values_to = "log2_vs_own_sex_vehicle") |>
    mutate(detected = detected[cbind(feature, sample_id)]) |>
    left_join(select(selected, feature, gene, direction, cluster), by = "feature") |>
    left_join(select(targets, sample_id, group), by = "sample_id") |>
    relocate(gene, direction, cluster, .after = feature) |>
    relocate(group, .after = sample_id)
  list(
    plot = wrap_elements(full = column) + theme_figure(),
    data = list(
      volcano = transmute(proteins,
        feature, gene, logFC,
        p = P.Value, fdr = adj.P.Val, pi_score, status
      ),
      heatmap = cells,
      ora = transmute(top,
        direction = selected$direction[match(cluster, selected$cluster)], cluster, collection,
        set = term, label, overlap, size, p, fdr
      )
    ),
    width_mm = width_mm, height_mm = height_mm
  )
}

# The key S3's columns share, in one row as a grob: the volcano points, the heatmap fill and
# not-detected grey, the direction strip, and the databases and FDR outline of the ORA bars.
pi_heatmap_key <- function(collections) {
  small <- grid::gpar(fontsize = 5.5)
  bold <- grid::gpar(fontsize = 5.5, fontface = "bold")
  collections <- intersect(names(DATABASE_COLOURS), collections)
  keys <- list(
    ComplexHeatmap::Legend(
      labels = VOLCANO_STATUS, type = "points", pch = 16, background = NA,
      legend_gp = grid::gpar(col = VOLCANO_COLOURS), title = "Volcano", title_gp = bold,
      labels_gp = small, nrow = 1, grid_width = unit(2, "mm"), grid_height = unit(2, "mm")
    ),
    ComplexHeatmap::Legend(
      col_fun = heatmap_colours, at = c(-HEATMAP_LIMIT, 0, HEATMAP_LIMIT),
      title = "log2 vs own-sex VEH", title_gp = bold, labels_gp = small,
      grid_height = unit(2, "mm"), legend_width = unit(14, "mm"), direction = "horizontal"
    ),
    ComplexHeatmap::Legend(
      labels = "not detected in that mouse", legend_gp = grid::gpar(fill = "grey75"),
      title = "Grey cells", title_gp = bold, labels_gp = small,
      grid_width = unit(2, "mm"), grid_height = unit(2, "mm")
    ),
    ComplexHeatmap::Legend(
      labels = c("up, or F > M", "down, or F < M"),
      legend_gp = grid::gpar(fill = DIRECTION_COLOURS[c("Up", "Down")]), title = "Direction",
      title_gp = bold, labels_gp = small, nrow = 1,
      grid_width = unit(2, "mm"), grid_height = unit(2, "mm")
    ),
    ComplexHeatmap::Legend(
      labels = collections, legend_gp = grid::gpar(fill = DATABASE_COLOURS[collections]),
      title = "ORA bars (-log10 p; dotted: p = 0.05)", title_gp = bold, labels_gp = small, nrow = 1,
      grid_width = unit(2, "mm"), grid_height = unit(2, "mm")
    ),
    ComplexHeatmap::Legend(
      labels = paste("FDR ≤", ALPHA), legend_gp = grid::gpar(fill = "white"), border = "black",
      title = "Outlined", title_gp = bold, labels_gp = small,
      grid_width = unit(2, "mm"), grid_height = unit(2, "mm")
    )
  )
  pack <- ComplexHeatmap::packLegend(
    list = keys, direction = "horizontal", column_gap = unit(5, "mm")
  )
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  grid::grid.grabExpr(ComplexHeatmap::draw(pack))
}
