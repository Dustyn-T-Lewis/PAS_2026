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
# Up and down in green and magenta (PiYG), so direction never borrows a group's orange or blue;
# fold changes, NES, correlations and r-equivalents all use it.
DIRECTION_COLOURS <- c(Down = "#C51B7D", Up = "#4D9221")
DIRECTION_LIGHT <- c(Down = "#F1B6DA", Up = "#B8E186")
PHENOTYPE_LABELS <- c(
  csa_mean = "CSA", capillary_density = "Capillary density",
  satellite_cells_per_fiber = "Satellite cells / fibre",
  conductance_sol_pyr = "Soleus pyruvate", conductance_pla_pyr = "Plantaris pyruvate",
  conductance_sol_oct = "Soleus octanoyl"
)
# Short names, for dense matrix columns only.
TRAIT_SHORT <- c(
  csa_mean = "CSA", capillary_density = "Capillaries", satellite_cells_per_fiber = "Sat. cells",
  conductance_sol_pyr = "Sol. pyr", conductance_pla_pyr = "Pla. pyr",
  conductance_sol_oct = "Sol. oct"
)
# One colour per phenotype for loading arrows, away from the group oranges and blues.
PHENOTYPE_COLOURS <- set_names(
  c("#1B9E77", "#7570B3", "#E6AB02", "#A6761D", "#666666", "#111111"), PHENOTYPE_LABELS
)
MODEL_LABELS <- c(pooled = "Pooled", treatment = "Within treatment")
# Within treatment gets a hue of its own, away from the interaction contrast's purple.
MODEL_COLOURS <- c(pooled = "grey40", pooled_mito = "grey70", treatment = "#B2182B")
CONTRAST_LABELS <- c(
  PAS_vs_VEH_female = "PAS vs VEH, female", PAS_vs_VEH_male = "PAS vs VEH, male",
  treatment_by_sex = "Treatment × sex", sex_in_VEH = "Sex in VEH", sex_in_PAS = "Sex in PAS"
)
CONTRAST_SHORT <- c(
  PAS_vs_VEH_female = "Female", PAS_vs_VEH_male = "Male", treatment_by_sex = "Treatment × sex"
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
# Paul Tol's muted and medium-contrast tones, away from the up and down magenta and green and the
# female orange and male blue, so a database never reads as a direction or a sex.
DATABASE_COLOURS <- c(
  "GO:BP" = "#44AA99", Reactome = "#332288", Hallmark = "#DDCC77", KEGG = "#88CCEE",
  MitoCarta = "#997700"
)
# WGCNA's module names, each drawn in a muted tone of its colour.
MODULE_COLOURS <- c(
  turquoise = "#8CCFC9", blue = "#8FA8D8", brown = "#BC9A7C", yellow = "#E6D58C",
  green = "#9DCB92", red = "#E39B94", black = "#8A8A8A", pink = "#E9B7CF"
)
# The response patterns of the moderated F over both PAS contrasts, in display order.
FTEST_PATTERNS <- c("up in both", "down in both", "female up, male down", "female down, male up")
# Dot area stops growing at this -log10 FDR or p.
LOG_CAP <- 4
# The every-mouse panels show the candidates passing at least this many routes; the other
# candidate thresholds come with candidates.rds from 02_Statistics/06_Candidates.
CANDIDATE_DOT_ROUTES <- 3

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
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.key.size = unit(2.5, "mm"),
      legend.text = element_text(size = 6, margin = margin(l = 1, r = 4)),
      # Room between a key's title and its labels (a colour bar's end label sits under its edge),
      # and between keys stacked in one box.
      legend.title = element_text(margin = margin(r = 8, b = 4)),
      legend.margin = margin(0), legend.box.spacing = unit(1.5, "mm"),
      legend.spacing.y = unit(2.5, "mm"),
      plot.tag = element_text(face = "bold", size = 9, vjust = 1, margin = margin(r = 2)),
      plot.tag.location = "margin", plot.tag.position = "left"
    )
}

# The key every continuous fill or colour scale uses.
colour_bar <- \(...) {
  guide_colourbar(
    barwidth = unit(14, "mm"), barheight = unit(2, "mm"),
    theme = theme(legend.text = element_text(margin = margin(t = 1.5))), ...
  )
}

# p as the labels print it, without scientific notation.
format_p <- \(p) if_else(p < 0.001, "p < 0.001", sprintf("p = %.2g", p))

# Gene symbols, with the accession added where two protein groups share one.
protein_label <- function(gene, feature) {
  shared <- gene %in% gene[duplicated(gene)]
  if_else(shared, paste(gene, feature), gene)
}

# Set labels, with the database added where two sets from different databases share one.
label_with_database <- function(label, collection) {
  shared <- label %in% label[duplicated(label)]
  if_else(shared, sprintf("%s (%s)", label, collection), label)
}

# Facet strips as the design panel's boxes: a tint inside a border in its colour, the name in bold
# black. Each strip is keyed, in display order, by a contrast or "MitoCarta" (its database colour);
# "" draws no box and any other key a grey one.
contrast_strips <- function(x = NULL, y = NULL) {
  fill <- c(contrast_tone("fill"), MitoCarta = "grey90")
  line <- c(contrast_tone("line"), DATABASE_COLOURS["MitoCarta"])
  boxes <- \(keys) {
    tone <- \(tones, other) {
      unname(if_else(keys %in% "", "transparent", coalesce(tones[keys], other)))
    }
    ggh4x::elem_list_rect(
      fill = tone(fill, "grey94"), colour = tone(line, "grey60"), linewidth = 0.4
    )
  }
  text <- ggh4x::elem_list_text(colour = "black", face = "bold")
  ggh4x::strip_themed(
    background_x = if (length(x)) boxes(x), text_x = if (length(x)) text,
    background_y = if (length(y)) boxes(y), text_y = if (length(y)) text
  )
}

# A figure panel at the size it takes in its composite, as PDF (cairo) and PNG, into the figure's
# b_reports/panels/; the composite then stitches panels that already look final.
save_panel <- function(plot, figure, name, width_mm, height_mm) {
  dir <- here("04_Figures", figure, "b_reports", "panels")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  save_figure(plot, dir, name, width_mm, height_mm)
}

# Mice as points, fill by group and shape by treatment (circle vehicle, triangle PAS), with one
# key that shows both; the points need shape = treatment, "VEH" or "PAS".
group_point_scales <- function() {
  list(
    scale_fill_manual(
      values = GROUP_COLOURS, labels = GROUP_LABELS, name = NULL,
      guide = guide_legend(
        nrow = 2, byrow = TRUE, order = 1, override.aes = list(shape = c(21, 24, 21, 24))
      )
    ),
    scale_shape_manual(values = c(VEH = 21, PAS = 24), guide = "none")
  )
}

# PCA scores with a shaded 95% t-ellipse per group (fill: group; shape: treatment), the PERMANOVA
# terms in the top-right corner (bold at p <= 0.05), and arrows (x, y, label) from the origin,
# already scaled to the score axes, each in its own colour with a key.
pca_panel <- function(scores, variance, title, subtitle, permanova, arrows, arrow_colours) {
  axis <- \(i) sprintf("PC%d (%.1f%%)", i, 100 * variance[i])
  scores <- mutate(scores,
    group = factor(group, FIGURE_GROUPS), treatment = factor(treatment, c("VEH", "PAS"))
  )
  arrows <- mutate(arrows, label = factor(label, names(arrow_colours)))
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
    geom_point(aes(fill = group, shape = treatment), size = 1.5, colour = "grey15", stroke = 0.25) +
    geom_text(
      aes(Inf, Inf, label = text, fontface = face, vjust = vjust), permanova,
      hjust = 1.04, size = 1.6, colour = "grey15", inherit.aes = FALSE
    ) +
    group_point_scales() +
    scale_colour_manual(
      values = arrow_colours, labels = scales::label_wrap(20), name = NULL,
      guide = guide_legend(ncol = 2, order = 2, override.aes = list(linewidth = 0.5))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    labs(x = axis(1), y = axis(2), title = title, subtitle = subtitle) +
    theme_figure() +
    theme(legend.box = "vertical", legend.key.width = unit(4, "mm"))
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

# Dot area grows with -log10 FDR (or p), capped at cap, from a minimum big enough to show the
# fill of a cell near FDR 1; the key reads in FDR (or p) units.
log_size_scale <- function(name = "FDR", range = c(0.6, 3.2), cap = LOG_CAP) {
  scale_size(
    range = range, limits = c(0, cap), breaks = -log10(c(0.1, 0.05, 0.01, 10^-cap)),
    labels = c("0.1", "0.05", "0.01", paste("≤", format(10^-cap, scientific = FALSE))),
    name = name
  )
}

# Dot size and ring shared by the dot matrices: log_size_scale(), and a thick ring at FDR <= ALPHA.
fdr_scales <- function(cap = LOG_CAP, name = "FDR") {
  list(
    log_size_scale(name, cap = cap),
    scale_discrete_manual(
      "stroke",
      values = c(`TRUE` = 0.7, `FALSE` = 0.15), breaks = "TRUE", labels = paste("FDR ≤", ALPHA),
      name = NULL, guide = guide_legend(override.aes = list(size = 2.4, fill = "white"))
    )
  )
}

# Pathway sets as a dot grid, after HRvLR's pathway_grid(): one row per set within its block, one
# column per contrast or trait under boxed strips, a dot per cell (fill value, keyed as fill_name
# on a PiYG scale to +-fill_limit, by default the largest drawn; area -log10 p, keyed as size_name)
# and, when rows has a collection, the set's database as a square before the first column. rows
# holds block (a factor in display order), term and label, top to bottom, and optionally
# collection and shade (a colour for a band behind the row's dots, NA for none); cells holds term,
# column, value and p (an FDR or a p), and optionally selected (TRUE outlines the cell in
# outline_shape, keyed as selected_label); columns holds column, column_block and label, in display
# order, the column blocks side by side under one key. tones maps a strip's name to its
# contrast_strips() key. Panels have a fixed size, column_mm wide and 3 mm a row.
pathway_grid <- function(cells, rows, columns, title = NULL, subtitle = NULL, tones = NULL,
                         x_angle = 0, column_mm = 5, selected_label = "Selected",
                         fill_name = "NES", fill_limit = NULL, fill_breaks = c(-2, 0, 2),
                         size_name = "FDR", outline_shape = 22) {
  column_blocks <- unique(columns$column_block)
  blocks <- levels(droplevels(rows$block))
  # The database squares take a narrow column of their own, " " under an empty strip.
  lead <- if (has_name(rows, "collection")) " " else character()
  place <- \(d) {
    mutate(d,
      term = factor(term, rev(rows$term)),
      column = factor(column, c(lead, columns$column)),
      column_block = factor(column_block, c(if (length(lead)) "", column_blocks))
    )
  }
  dots <- cells |>
    inner_join(select(rows, term, block, any_of("shade")), by = "term") |>
    inner_join(select(columns, column, column_block), by = "column") |>
    mutate(log_p = pmin(-log10(p), LOG_CAP)) |>
    place()
  # One full-width tile per row and panel, so the alpha does not stack across columns.
  shaded <- dots |>
    filter(!is.na(dots$shade %||% NA)) |>
    distinct(term, column_block, .keep_all = TRUE)
  band <- if (nrow(shaded)) {
    geom_tile(data = shaded, fill = shaded$shade, width = Inf, height = 0.86)
  }
  outlined <- filter(dots, dots$selected %||% FALSE)
  squares <- if (length(lead)) {
    geom_point(
      aes(colour = collection), place(mutate(rows, column = lead, column_block = "")),
      shape = 15, size = 1.8, show.legend = c(colour = TRUE, size = FALSE, fill = FALSE)
    )
  }
  fill_limit <- fill_limit %||% (ceiling(max(abs(dots$value)) * 10) / 10)
  ggplot(dots, aes(column, term)) +
    band +
    geom_point(
      aes(shape = selected_label), outlined,
      size = 3.6, fill = NA, colour = "grey15", stroke = 0.35
    ) +
    geom_point(aes(size = log_p, fill = value), shape = 21, stroke = 0.2) +
    squares +
    ggh4x::facet_grid2(block ~ column_block,
      scales = "free", space = "free", switch = "y",
      strip = contrast_strips(c(if (length(lead)) "", tones[column_blocks]), tones[blocks])
    ) +
    ggh4x::force_panelsizes(
      cols = unit(
        c(if (length(lead)) 3.5, column_mm * table(columns$column_block)[column_blocks]), "mm"
      ),
      rows = unit(3 * table(rows$block)[blocks], "mm")
    ) +
    # Order pinned from rows and columns, so a layer that lacks some rows (the bands, the
    # outlines) cannot reorder a facet and move an outline off its cell.
    scale_y_discrete(
      limits = \(x) intersect(rev(rows$term), x), labels = deframe(select(rows, term, label)),
      expand = expansion(add = 0.6)
    ) +
    scale_x_discrete(
      limits = \(x) intersect(c(lead, columns$column), x),
      labels = deframe(select(columns, column, label))
    ) +
    scale_fill_distiller(
      palette = "PiYG", direction = 1, limits = c(-fill_limit, fill_limit), oob = scales::squish,
      breaks = fill_breaks, name = fill_name
    ) +
    log_size_scale(size_name, range = c(0.6, 2.8)) +
    scale_colour_manual(
      values = DATABASE_COLOURS, limits = names(DATABASE_COLOURS), name = "Database"
    ) +
    scale_shape_manual(values = outline_shape, name = NULL) +
    guides(
      fill = guide_colourbar(order = 1, barheight = unit(12, "mm"), barwidth = unit(2.2, "mm")),
      size = guide_legend(order = 2), colour = guide_legend(order = 3),
      shape = guide_legend(order = 4, override.aes = list(size = 3))
    ) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_figure() +
    theme(
      axis.text.y = element_text(size = 5.2),
      axis.text.x = element_text(
        size = 5.5, angle = x_angle, hjust = if (x_angle) 1 else 0.5, vjust = 0.5
      ),
      strip.placement = "outside", strip.text = element_text(size = 5.5),
      strip.text.y.left = element_text(angle = 90), panel.spacing = unit(1, "mm"),
      legend.position = "right", legend.box.just = "left", legend.key.size = unit(2.6, "mm"),
      legend.text = element_text(size = 5.5), legend.title = element_text(size = 6)
    )
}

# Column strips of the F03 and F04 grids, keyed to contrast_strips(); Contrast falls back to grey.
SEX_TRAIT_TONES <- c(`Female traits` = "PAS_vs_VEH_female", `Male traits` = "PAS_vs_VEH_male")

# The columns F03 and F04 share, for the sets in terms: the fgsea result in each contrast, then
# the within-treatment association fgsea at each trait in females and in males, one row per cell
# in column order. column, column_block and label are pathway_grid()'s columns; fdr is BH.
contrast_trait_cells <- function(enr, terms) {
  contrasts <- enr$results |>
    filter(term %in% terms) |>
    transmute(
      column_block = "Contrast", column = contrast, label = CONTRAST_SHORT[contrast], contrast,
      term, collection, universe, NES, p, fdr = padj, shown, camera_p, camera_fdr
    )
  traits <- enr$association_fgsea |>
    filter(model == "treatment", term %in% terms) |>
    transmute(
      column_block = names(SEX_TRAIT_TONES)[match(sex, c("female", "male"))],
      column = paste(sex, phenotype), label = TRAIT_SHORT[phenotype], sex, phenotype, term,
      collection, universe, NES, p = fgsea_p, fdr = fgsea_fdr, shown, camera_p, camera_fdr
    )
  bind_rows(contrasts, traits) |>
    arrange(
      factor(column_block, c("Contrast", names(SEX_TRAIT_TONES))),
      match(contrast, COMPARISON_CONTRASTS), match(phenotype, names(TRAIT_SHORT))
    ) |>
    relocate(sex, phenotype, .after = contrast)
}

# Each set once per figure: block by block in factor order, the first n candidates (already in
# display order) that no earlier block lists.
list_once <- function(candidates, n) {
  reduce(levels(candidates$block), \(listed, b) {
    new <- candidates |>
      filter(block == b, !term %in% listed$term) |>
      slice_head(n = n)
    bind_rows(listed, new)
  }, .init = candidates[0, ])
}

# Set counts beside a pathway_grid(), on its row bands and blocks: per group, an up bar above a
# down bar from zero, each stacked by database inside one black outline, its total at the end, a
# ▲ or ▼ before it and the group's label before the pair. pairs holds block, bar (the label), tone
# (its colour), y (the pair's centre in rows from the block's top) and span (the rows it covers);
# counts holds block, bar, direction ("up" or "down"), collection and n; rows is the grid's rows.
count_bars <- function(counts, pairs, rows, width_mm = 30,
                       x_title = sprintf(
                         "Distinct significant sets\n(FDR ≤ %s, after overlap grouping)", ALPHA
                       )) {
  blocks <- levels(droplevels(rows$block))
  sizes <- table(rows$block)[blocks]
  gap <- 0.15
  # Rows flipped onto a continuous axis, so a centre at k sits level with the grid's k-th row; a
  # bar is up to 1.4 rows thick and a pair fills at most 90% of its rows.
  pairs <- mutate(pairs,
    block = factor(block, blocks), y = as.numeric(sizes[block]) + 1 - y,
    thick = pmin(1.4, (0.9 * span - gap) / 2)
  )
  slots <- pairs |>
    cross_join(tibble(direction = c("down", "up"))) |>
    mutate(centre = y + if_else(direction == "up", 1, -1) * (thick + gap) / 2)
  bars <- counts |>
    filter(n > 0) |>
    mutate(
      block = factor(block, blocks), collection = factor(collection, names(DATABASE_COLOURS))
    ) |>
    inner_join(
      select(slots, block, bar, direction, centre, thick),
      by = c("block", "bar", "direction")
    ) |>
    arrange(collection) |>
    mutate(xmax = cumsum(n), xmin = xmax - n, .by = c(block, bar, direction))
  # Every slot gets its glyph and total, 0 where the direction has no sets.
  ends <- slots |>
    left_join(
      summarise(bars, x = sum(n), .by = c(block, bar, direction)),
      by = c("block", "bar", "direction")
    ) |>
    mutate(x = coalesce(x, 0L))
  reach <- max(ends$x) * 1.2
  mm <- reach / width_mm
  bands <- rows |>
    mutate(y = n() + 1 - row_number(), .by = block) |>
    filter(!is.na(shade))
  # The grid's 0.6-row padding above and below each block.
  frame <- tibble(
    block = factor(rep(blocks, 2), blocks), y = c(rep(0.4, length(blocks)), sizes + 0.6)
  )
  ggplot(bars) +
    geom_blank(aes(0, y), frame) +
    geom_tile(aes(0, y), bands, fill = bands$shade, width = Inf, height = 0.86) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey40") +
    geom_rect(aes(
      xmin = xmin, xmax = xmax, ymin = centre - thick / 2, ymax = centre + thick / 2,
      fill = collection
    )) +
    # One outline per bar, not per slice, so a one-set slice keeps its colour.
    geom_rect(
      aes(xmin = 0, xmax = x, ymin = centre - thick / 2, ymax = centre + thick / 2),
      filter(ends, x > 0),
      fill = NA, colour = "black", linewidth = 0.33
    ) +
    geom_text(aes(x + mm, centre, label = x), ends, hjust = 0, size = 1.9) +
    geom_text(
      aes(-1.6 * mm, centre,
        label = if_else(direction == "up", "▲", "▼"),
        colour = DIRECTION_COLOURS[if_else(direction == "up", "Up", "Down")]
      ), ends,
      size = 1.5
    ) +
    geom_text(aes(-3.4 * mm, y, label = bar, colour = tone), pairs,
      hjust = 1, size = 1.7, fontface = "bold"
    ) +
    facet_grid(block ~ ., scales = "free_y", space = "free_y") +
    ggh4x::force_panelsizes(cols = unit(width_mm, "mm"), rows = unit(3 * sizes, "mm")) +
    scale_x_continuous(breaks = scales::breaks_extended(5)) +
    scale_fill_manual(values = DATABASE_COLOURS, guide = "none") +
    scale_colour_identity() +
    coord_cartesian(xlim = c(0, reach), expand = FALSE, clip = "off") +
    labs(x = x_title, y = NULL) +
    theme_figure() +
    theme(
      axis.text.y = element_blank(), panel.grid.major.y = element_blank(),
      strip.text = element_blank(), panel.spacing = unit(1, "mm"),
      axis.text.x = element_text(size = 5), axis.title.x = element_text(size = 5.5),
      plot.margin = margin(l = 22, r = 4)
    )
}

# Display names for pathway sets: the collection prefix dropped, title case, shortened.
set_label <- \(term, width = 48) str_trunc(enrichVolcano::clean_label(term, width = 200), width)

# protein_heatmap()'s fill, full magenta at -limit and full green at +limit, and its key: the fill
# with breaks at -limit, 0 and limit, then the not-detected grey.
heatmap_colours <- \(limit) {
  circlize::colorRamp2(
    c(-limit, 0, limit), c(DIRECTION_COLOURS[["Down"]], "white", DIRECTION_COLOURS[["Up"]])
  )
}
heatmap_legends <- function(limit, name = "log2 vs\nown-sex VEH", direction = "vertical") {
  small <- grid::gpar(fontsize = 5.5)
  list(
    ComplexHeatmap::Legend(
      col_fun = heatmap_colours(limit), at = c(-limit, 0, limit), title = name,
      title_gp = grid::gpar(fontsize = 5.5, fontface = "bold"), labels_gp = small,
      grid_width = unit(2, "mm"), grid_height = unit(2, "mm"), legend_width = unit(14, "mm"),
      direction = direction
    ),
    ComplexHeatmap::Legend(
      labels = "not detected", legend_gp = grid::gpar(fill = "grey75"), labels_gp = small,
      grid_width = unit(2, "mm"), grid_height = unit(2, "mm")
    )
  )
}

# Row order for a protein heatmap: coexpression_tree() (helpers.R) within each block (blocks, a
# factor, in level order) of the drawn values x.
coexpression_order <- function(x, blocks, min_shared = 5) {
  map(split(seq_len(nrow(x)), blocks), \(rows) {
    rows[coexpression_tree(x[rows, , drop = FALSE], min_shared)$order]
  }) |>
    unlist(use.names = FALSE)
}

# Proteins by mice as a ComplexHeatmap, each protein centred on its own sex's vehicle mean, rows
# drawn in the order given (see coexpression_order()), columns split by group under a group-colour
# bar and clustered within it, magenta below and green above; cells where detected (same shape as
# values) is FALSE are grey, so an imputed value never reads as no change. targets lists the
# columns of values in order. Each row_split level is labelled with its protein count after
# split_sep, unless row_titles is FALSE; row_gap (mm, one per gap or one for all) separates the
# levels. Returns a grob patchwork can place with wrap_elements(); left and right are optional
# rowAnnotations; padding (mm: bottom, left, top, right) surrounds the heatmap. Group titles take
# two lines, or one line upright with the block titles when rotate_titles is TRUE, for a narrow
# heatmap.
protein_heatmap <- function(values, targets, detected, row_split = NULL, split_sep = " ",
                            row_titles = TRUE, row_gap = 0.8, left = NULL, right = NULL,
                            limit = 1, name = "log2 vs\nown-sex VEH", width_mm = 120,
                            height_mm = 120, legend = TRUE, rotate_titles = FALSE,
                            padding = c(1, 1, 4, 1)) {
  groups <- factor(targets$group, FIGURE_GROUPS)
  group_titles <- GROUP_LABELS[levels(droplevels(groups))]
  if (!is.null(row_split) && row_titles) {
    row_split <- droplevels(row_split)
    row_split <- factor(row_split,
      labels = sprintf("%s%s(%d)", levels(row_split), split_sep, tabulate(row_split))
    )
  }
  top <- ComplexHeatmap::HeatmapAnnotation(
    group = groups, col = list(group = GROUP_COLOURS), show_annotation_name = FALSE,
    show_legend = FALSE, simple_anno_size = unit(1.6, "mm")
  )
  small <- grid::gpar(fontsize = 5.5)
  not_detected <- \(j, i, x, y, w, h, fill) {
    hide <- !ComplexHeatmap::pindex(detected, i, j)
    if (any(hide)) {
      grid::grid.rect(x[hide], y[hide], w[hide], h[hide],
        gp = grid::gpar(fill = "grey75", col = NA)
      )
    }
  }
  heatmap <- ComplexHeatmap::Heatmap(centre_on_vehicle(values, targets),
    name = name, col = heatmap_colours(limit), top_annotation = top, left_annotation = left,
    right_annotation = right, column_split = groups, cluster_column_slices = FALSE,
    column_title_gp = small,
    column_title = if (rotate_titles) group_titles else str_replace(group_titles, " ", "\n"),
    column_title_rot = if (rotate_titles) 90 else 0, column_labels = rep("", ncol(values)),
    row_split = row_split, cluster_rows = FALSE, row_title_gp = small,
    row_title = if (row_titles) character(0),
    row_title_rot = if (rotate_titles) 90 else 0, show_row_names = FALSE,
    show_column_dend = FALSE, column_gap = unit(0.6, "mm"), row_gap = unit(row_gap, "mm"),
    layer_fun = not_detected, show_heatmap_legend = FALSE
  )
  # grid.grabExpr returns to the device open before it, and with none open R starts Rplots.pdf.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  grid::grid.grabExpr(
    ComplexHeatmap::draw(heatmap,
      heatmap_legend_side = "right", annotation_legend_side = "right",
      heatmap_legend_list = if (legend) heatmap_legends(limit, name),
      padding = unit(padding, "mm")
    ),
    width = width_mm / 25.4, height = height_mm / 25.4
  )
}

# A figure's panel scripts, each run in its own environment, as YvO's stitchers do, and each plot
# given its letter when there is more than one. A panel script returns list(plot, data, inputs).
run_panels <- function(figure, panels) {
  map(set_names(panels), \(name) {
    path <- here("04_Figures", figure, "a_script", "panels", paste0(name, ".R"))
    panel <- source(path, local = new.env())$value
    if (length(panels) > 1) panel$plot <- panel$plot + labs(tag = letters[match(name, panels)])
    panel
  })
}

# A composite's workbook: one sheet per panel (a panel whose data is a list of tables gives one
# sheet per table, named <panel>_<table>), then every input the panels read.
write_figure_workbook <- function(panels, holds, file) {
  inputs <- unlist(unname(map(panels, "inputs")))
  sheets <- list_flatten(map(panels, "data"), name_spec = "{outer}_{inner}")
  write_workbook(sheets, holds, file, inputs[!duplicated(inputs)])
}

# One contrast's S2 heatmap: its proteins at Π <= 0.05 (full proteome) across every mouse, each
# centred on its own sex's vehicle mean, grey where not detected, in 02_Enrichment's drawing order
# (heatmap_clusters): split by direction (for the interaction, by which sex's PAS effect is
# larger) and within each direction by co-expression cluster, the proteins the cut left
# unassigned at the block's end under its last cluster. Beside each cluster, up to n_sets sets its
# proteins over-represent at p <= ALPHA, smallest p (then FDR) first (cluster_ora), on one -log10 p
# range for every comparison contrast. fit is 01_Fit's output, quant 02_Quantification's and
# enrichment 02_Enrichment's. width_mm and height_mm are the panel's size, returned for
# save_panel; the title and subtitle take the top 10 mm and the bars ora_mm of the width.
pi_heatmap_panel <- function(fit, quant, enrichment, contrast, width_mm = 60, height_mm = 178,
                             ora_mm = 22, n_sets = 3) {
  names <- if (contrast == "treatment_by_sex") c("F > M", "F < M") else c("up", "down")
  targets <- arrange(fit$targets, factor(group, FIGURE_GROUPS))
  selected <- enrichment$heatmap_clusters |>
    filter(contrast == !!contrast) |>
    arrange(row) |>
    left_join(
      filter(fit$protein_contrasts, universe == "full", contrast == !!contrast) |>
        select(feature, logFC, P.Value, adj.P.Val, pi_score),
      by = "feature"
    ) |>
    mutate(direction = factor(names[match(direction, c("up", "down"))], names), slice = cluster) |>
    fill(slice, .by = direction)
  slices <- factor(selected$slice, unique(selected$slice))
  values <- fit$E_norm$full[selected$feature, targets$sample_id]
  detected <- quant$proteins$other$n.observations[selected$feature, targets$sample_id] > 0
  top <- enrichment$cluster_ora |>
    filter(p <= ALPHA) |>
    arrange(p, fdr) |>
    slice_head(n = n_sets, by = c(contrast, cluster))
  reach <- max(-log10(c(top$p, ALPHA))) * 1.05
  top <- top |>
    filter(contrast == !!contrast) |>
    arrange(factor(cluster, levels(slices)))
  block_of <- map_chr(split(as.character(selected$direction), slices), 1)
  gaps <- if_else(head(block_of, -1) == tail(block_of, -1), 0.7, 2.6)
  # A direction block's strip, drawn from its first slice down across its slices and gaps.
  strip <- function(index, levels) {
    first <- as.integer(slices[index[1]])
    if (first > 1 && block_of[[first - 1]] == block_of[[first]]) {
      return(invisible())
    }
    at <- slice_positions(slices, gaps, length(index))
    block <- block_of[[first]]
    span <- max((at$top + at$height)[block_of == block]) - at$top[first]
    # In npc of this slice, so the strip scales with the grob when patchwork fits it to the panel.
    grid::pushViewport(grid::viewport(
      y = unit(1, "npc"), height = unit(span / at$height[first], "npc"), just = "top"
    ))
    tone <- if (block == names[1]) "Up" else "Down"
    grid::grid.rect(gp = grid::gpar(fill = DIRECTION_LIGHT[[tone]], col = NA))
    grid::grid.text(
      sprintf("%s (%d)", block, sum(selected$direction == block)),
      rot = 90, gp = grid::gpar(fontsize = 5.5)
    )
    grid::popViewport()
  }
  cluster_n <- table(factor(selected$cluster, levels(slices)))
  left <- ComplexHeatmap::rowAnnotation(
    direction = ComplexHeatmap::anno_block(
      which = "row", width = unit(3, "mm"), panel_fun = strip
    ),
    cluster = ComplexHeatmap::anno_block(
      which = "row", width = unit(7.5, "mm"), gp = grid::gpar(fill = NA, col = NA),
      labels = sprintf("%s (%d)", names(cluster_n), cluster_n), labels_rot = 0,
      labels_just = "right", labels_offset = unit(1, "npc") - unit(0.6, "mm"),
      labels_gp = grid::gpar(fontsize = 5)
    ),
    show_annotation_name = FALSE, gap = unit(0.4, "mm")
  )
  right <- ComplexHeatmap::rowAnnotation(
    ora = ora_bars(top, slices, block_of, gaps, ora_mm, reach), show_annotation_name = FALSE
  )
  heatmap <- protein_heatmap(values, targets, detected,
    row_split = slices, left = left, right = right, row_titles = FALSE, row_gap = gaps,
    limit = 1.5, width_mm = width_mm, height_mm = height_mm - 10, legend = FALSE,
    rotate_titles = TRUE, padding = c(7, 1, 4, 1)
  )
  plot <- wrap_elements(plot = heatmap) +
    labs(
      title = CONTRAST_LABELS[[contrast]],
      subtitle = sprintf("Π ≤ %s: %d proteins; bars: ORA", ALPHA, nrow(selected))
    ) +
    theme_figure()
  list(
    plot = plot,
    data = list(
      heatmap = select(
        selected, feature, gene, direction, cluster, logFC, P.Value, adj.P.Val, pi_score
      ),
      ora = select(
        top, contrast, direction, cluster, collection, term, overlap, size, p, fdr,
        genes_in_cluster, overlap_genes
      )
    ),
    width_mm = width_mm, height_mm = height_mm
  )
}

# Set names for the bars beside a heatmap: sentence case (proper nouns kept), a few long words
# shortened, and the hyphens clean_label() drops put back.
PROPER_WORDS <- "Golgi"
WORD_ABBREVIATIONS <- c(
  regulation = "reg.", development = "dev.", process = "proc.", organization = "org.",
  positive = "pos.", negative = "neg.", dependent = "dep.", localization = "loc.",
  biosynthetic = "biosynth."
)
NAME_FIXES <- c("\\b([Cc])ell cell\\b" = "\\1ell-cell", " (mediated|based|dep\\.)(?= |$)" = "-\\1")
DATABASE_TAGS <- c("GO:BP" = "GO", Reactome = "R", Hallmark = "H", KEGG = "K")
bar_set_name <- function(term) {
  words <- str_split(set_label(term, 200), " ")
  map_chr(words, \(w) {
    lower <- seq_along(w) > 1 & str_detect(w, "^[A-Z][a-z]+$") & !w %in% PROPER_WORDS
    w[lower] <- str_to_lower(w[lower])
    short <- WORD_ABBREVIATIONS[str_to_lower(w)]
    w[!is.na(short)] <- paste0(str_sub(w, 1, 1), str_sub(short, 2))[!is.na(short)]
    paste(w, collapse = " ")
  }) |>
    str_replace_all(NAME_FIXES)
}

# Each row slice's top and height in mm down a split heatmap's body, from inside one slice's
# annotation viewport: that slice's rows (n_rows) fix the mm per row, gaps (mm) separate slices.
slice_positions <- function(slices, gaps, n_rows) {
  mm_per_row <- grid::convertHeight(unit(1, "npc"), "mm", valueOnly = TRUE) / n_rows
  height <- tabulate(slices) * mm_per_row
  list(top = cumsum(c(0, head(height, -1) + gaps)), height = height)
}

# Frames of the given sizes (mm) as near their centres (mm down the body) as they go without
# overlapping, none ending below its bottom (mm, one per frame): pushed down from the top, then
# back up from the bottom.
place_frames <- function(centre, size, bottom, pad = 2.5) {
  top <- centre - size / 2
  for (k in seq_along(top)[-1]) top[k] <- max(top[k], top[k - 1] + size[k - 1] + pad)
  for (k in rev(seq_along(top))) {
    limit <- min(bottom[k], if (k < length(top)) top[k + 1] - pad else Inf)
    top[k] <- min(top[k], limit - size[k])
  }
  top + size / 2
}

# The ORA bars beside each row slice of a heatmap, as a ComplexHeatmap row annotation: per slice,
# its sets stacked from the top in a frame beside the slice (centred, unless that would overlap a
# neighbouring frame or run past the slice's direction block), each a -log10 p bar bar_mm thick
# on the range 0 to reach, filled by database, outlined black at FDR <= ALPHA, with a dotted mark
# at p = ALPHA, and named "Set (overlap/size)" above it on at most two lines of `width` characters
# (a database tag where two of the slice's sets share a name; an ellipsis after the last whole
# word that fits when two lines cannot hold the name). A slice with no set says so. The axis
# hangs under the last slice. top holds cluster (the slice), collection, term, overlap, size, p
# and fdr, in drawing order; slices gives each heatmap row's slice, blocks each slice's direction
# block and gaps the heatmap's row gaps (mm).
ora_bars <- function(top, slices, blocks, gaps, width_mm, reach, bar_mm = 1.9, width = 24) {
  line_mm <- 1.75
  note_mm <- 2
  # The name's first k words, ellipsis after, for the largest k that fits two lines.
  two_lines <- function(name, counts) {
    words <- str_split_1(name, " ")
    for (k in rev(seq_along(words))) {
      name <- paste(head(words, k), collapse = " ")
      if (k < length(words)) name <- paste0(str_remove(name, "\\.$"), "…")
      label <- str_wrap(paste(name, counts), width)
      if (str_count(label, "\n") <= 1) break
    }
    label
  }
  top <- top |>
    mutate(
      name = bar_set_name(term),
      name = if_else(
        name %in% name[duplicated(name)], sprintf("%s (%s)", name, DATABASE_TAGS[collection]), name
      ),
      label = map2_chr(name, sprintf("(%d/%d)", overlap, size), two_lines),
      # Slot height: the label's lines, 0.8 mm to the bar, the bar and 0.8 mm below it.
      slot = (str_count(label, "\n") + 1) * line_mm + 0.8 + bar_mm + 0.8,
      bar_hi = sum(slot) - cumsum(slot) + 0.8 + bar_mm,
      .by = cluster
    )
  frame_mm <- map_dbl(levels(slices), \(s) max(sum(top$slot[top$cluster == s]), note_mm))
  bar_theme <- \() {
    theme_void(base_size = 7, base_family = "Helvetica") + theme(plot.margin = margin(0, 0, 0, 1))
  }
  bar_scales <- \(height) {
    list(
      scale_fill_manual(values = DATABASE_COLOURS, guide = "none"),
      scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "transparent"), guide = "none"),
      coord_cartesian(xlim = c(0, reach), ylim = c(0, height), expand = FALSE, clip = "off")
    )
  }
  # A name of its own per grob, so grid.grabExpr() keeps every slice's bars.
  draw_grob <- function(plot, y, height, just = "centre") {
    grid::pushViewport(grid::viewport(y = y, height = height, just = c("centre", just)))
    grob <- ggplotGrob(plot)
    grid::grid.draw(grid::editGrob(grob, name = grid::grobName(grob, "ora_bars")))
    grid::popViewport()
  }
  draw <- function(index, levels) {
    k <- as.integer(slices[index[1]])
    slice <- levels(slices)[k]
    at <- slice_positions(slices, gaps, length(index))
    centre <- at$top + at$height / 2
    # Each frame ends by the bottom of its direction block.
    block_bottom <- ave(at$top + at$height, blocks, FUN = max)
    placed <- place_frames(centre, frame_mm, block_bottom)
    y <- unit(0.5 - (placed[k] - centre[k]) / at$height[k], "npc")
    sets <- filter(top, cluster == slice)
    if (nrow(sets) > 0) {
      plot <- ggplot(sets) +
        geom_rect(aes(
          xmin = 0, xmax = -log10(p), ymin = bar_hi - bar_mm, ymax = bar_hi,
          fill = collection, colour = fdr <= ALPHA
        ), linewidth = 0.3) +
        geom_segment(aes(
          x = -log10(ALPHA), xend = -log10(ALPHA), y = bar_hi - bar_mm - 0.2, yend = bar_hi + 0.2
        ), linetype = "11", linewidth = 0.3) +
        geom_text(aes(0, bar_hi + 0.8, label = label),
          hjust = 0, vjust = 0, size = 5 / .pt, lineheight = 0.85, family = "Helvetica"
        ) +
        bar_scales(frame_mm[k]) +
        bar_theme()
      draw_grob(plot, y, unit(frame_mm[k], "mm"))
    } else {
      grid::grid.text(paste("no set at p ≤", ALPHA),
        x = unit(1, "mm"), y = y, just = "left",
        gp = grid::gpar(fontsize = 4.5, col = "grey50", fontfamily = "Helvetica")
      )
    }
    if (k == nlevels(slices)) {
      axis <- ggplot(tibble(x = c(0, reach)), aes(x, 1)) +
        geom_blank() +
        bar_scales(1) +
        scale_x_continuous(breaks = scales::breaks_pretty(4)) +
        labs(x = "-log10 p") +
        bar_theme() +
        theme(
          axis.line.x = element_line(linewidth = 0.25),
          axis.ticks.x = element_line(linewidth = 0.25),
          axis.ticks.length.x = unit(0.6, "mm"),
          axis.text.x = element_text(size = 5, margin = margin(t = 0.5)),
          axis.title.x = element_text(size = 5.5, margin = margin(t = 0.8))
        )
      draw_grob(axis, unit(0, "npc") - unit(1, "mm"), unit(6.5, "mm"), just = "top")
    }
  }
  ComplexHeatmap::anno_block(which = "row", width = unit(width_mm, "mm"), panel_fun = draw)
}

# The key S2's heatmaps share: the heatmap fill and not-detected grey, the databases its ORA bars
# use and the FDR outline, in one row, as a grob.
pi_heatmap_key <- function(collections, limit = 1.5) {
  small <- grid::gpar(fontsize = 5.5)
  collections <- intersect(names(DATABASE_COLOURS), collections)
  keys <- c(
    heatmap_legends(limit, "log2 vs own-sex VEH", direction = "horizontal"),
    list(
      ComplexHeatmap::Legend(
        labels = collections, legend_gp = grid::gpar(fill = DATABASE_COLOURS[collections]),
        title = "ORA database", title_gp = grid::gpar(fontsize = 5.5, fontface = "bold"),
        labels_gp = small, nrow = 1, grid_width = unit(2, "mm"), grid_height = unit(2, "mm")
      ),
      ComplexHeatmap::Legend(
        labels = paste("FDR ≤", ALPHA), legend_gp = grid::gpar(fill = "white"), border = "black",
        title = "Outlined", title_gp = grid::gpar(fontsize = 5.5, fontface = "bold"),
        labels_gp = small, grid_width = unit(2, "mm"), grid_height = unit(2, "mm")
      )
    )
  )
  pack <- ComplexHeatmap::packLegend(
    list = keys, direction = "horizontal", column_gap = unit(6, "mm")
  )
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  grid::grid.grabExpr(ComplexHeatmap::draw(pack))
}

# Sets by contrast as dots without a collection column, for side lists: NES as fill, -log10 FDR
# as size, rows grouped by row_group. cells holds row (display order), label, row_group, contrast,
# NES and fdr.
set_dot_matrix <- function(cells, title, subtitle, fdr_cap = LOG_CAP, nes_limit = 2.5) {
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
    scale_x_discrete(labels = scales::label_wrap(10)) +
    scale_y_discrete(labels = deframe(distinct(cells, row, label))) +
    scale_fill_gradient2(
      low = DIRECTION_COLOURS[["Down"]], high = DIRECTION_COLOURS[["Up"]],
      limits = c(-nes_limit, nes_limit), breaks = c(-2, 0, 2), name = "NES",
      guide = colour_bar()
    ) +
    fdr_scales(fdr_cap) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_figure() +
    theme(
      axis.text.y = element_text(size = 5.2), strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, size = 5.2),
      panel.spacing = unit(0.8, "mm"), legend.box = "vertical"
    )
}

# Each score's r-equivalent against the comparison contrasts and, within each sex, against each
# trait (pooled model, adjusted for the surrogate variable), as the cells of one matrix.
# contrasts and associations are 01_Fit's score_contrasts and score_associations, already
# filtered to the features shown.
score_cells <- function(contrasts, associations) {
  contrast_cells <- contrasts |>
    filter(contrast %in% COMPARISON_CONTRASTS) |>
    transmute(
      level, feature,
      section = "Contrast", column = CONTRAST_SHORT[contrast], r_equiv, P.Value, adj.P.Val
    )
  trait_cells <- associations |>
    filter(model == "pooled") |>
    transmute(
      level, feature,
      section = if_else(sex == "female", "Female traits", "Male traits"),
      column = TRAIT_SHORT[phenotype], r_equiv, pearson_r = r, P.Value, adj.P.Val
    )
  bind_rows(contrast_cells, trait_cells) |>
    mutate(
      section = factor(section, c("Contrast", "Female traits", "Male traits")),
      column = factor(column, c(CONTRAST_SHORT, TRAIT_SHORT))
    )
}

# One sex's candidates (those passing community_routes or more routes) as a STRING network:
# shaded hulls for communities, each named beside it by its ORA term; edges at the STRING score
# 06_Candidates used; node size routes passed, colour log2 fold change. Communities are kept
# together by adding invisible within-community edges to the layout only. Proteins without a
# STRING partner are left out of the drawing, counted under it and flagged by has_partner in the
# data. cand is candidates.rds from 02_Statistics/06_Candidates. Returns the plot and the
# candidates drawn.
candidate_network <- function(cand, sex, title) {
  limits <- cand$thresholds
  shown <- filter(cand$candidates, sex == !!sex, routes >= limits$community_routes)
  nodes <- shown |>
    distinct(gene, .keep_all = TRUE) |>
    mutate(group = if_else(is.na(name), NA_character_, str_wrap(set_label(name, 200), 24)))
  links <- cand$edges |>
    filter(sex == !!sex, a %in% nodes$gene, b %in% nodes$gene) |>
    distinct(a, b)
  linked <- nodes |> filter(!is.na(community) | gene %in% c(links$a, links$b))
  alone <- setdiff(nodes$gene, linked$gene)
  within <- linked |>
    filter(!is.na(community)) |>
    reframe(as_tibble(t(combn(gene, 2)), .name_repair = \(x) c("a", "b")), .by = community) |>
    select(a, b)
  layout_graph <- igraph::graph_from_data_frame(
    bind_rows(mutate(links, w = 1), mutate(within, w = 4)),
    directed = FALSE, vertices = select(linked, gene)
  )
  set.seed(SEED)
  xy <- igraph::layout_with_fr(layout_graph, weights = igraph::E(layout_graph)$w)
  positions <- tibble(gene = igraph::V(layout_graph)$name, x = xy[, 1], y = xy[, 2]) |>
    left_join(linked, by = "gene")
  segments <- links |>
    left_join(select(positions, a = gene, x, y), by = "a") |>
    left_join(select(positions, b = gene, xend = x, yend = y), by = "b")
  # Gene and community names repelled together, so neither covers a node or the other; each
  # community name starts above its hull's highest member.
  labels <- bind_rows(
    transmute(positions, x, y, label = gene, size = 1.9, colour = "black"),
    positions |>
      filter(!is.na(group)) |>
      summarise(x = mean(x), y = max(y), .by = group) |>
      transmute(x, y, label = group, size = 1.6, colour = "grey30")
  )
  limit <- max(abs(nodes$logFC))
  plot <- ggplot(positions, aes(x, y)) +
    ggforce::geom_mark_hull(
      aes(fill = group, filter = !is.na(group)),
      expand = unit(3.5, "mm"), radius = unit(3, "mm"), alpha = 0.18, colour = NA
    ) +
    geom_segment(aes(x, y, xend = xend, yend = yend), segments,
      colour = "grey60", linewidth = 0.3
    ) +
    geom_point(aes(colour = logFC, size = routes)) +
    ggrepel::geom_text_repel(aes(label = label), labels,
      size = labels$size, colour = labels$colour, lineheight = 0.9, seed = SEED,
      max.overlaps = Inf, box.padding = 0.15, point.padding = 0.1,
      min.segment.length = 0.3, segment.size = 0.2, bg.color = "white", bg.r = 0.15
    ) +
    scale_colour_gradient2(
      low = DIRECTION_COLOURS[["Down"]], high = DIRECTION_COLOURS[["Up"]],
      limits = c(-limit, limit), name = "log2 FC", guide = colour_bar(order = 1)
    ) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    # One size key for both sexes: the limits span every candidate passing community_routes.
    scale_size_continuous(
      range = c(1.6, 3.6), limits = c(limits$community_routes, max(cand$candidates$routes)),
      breaks = seq(limits$community_routes, max(cand$candidates$routes)), name = "routes",
      guide = guide_legend(order = 2)
    ) +
    # Room around the network for the community names.
    scale_x_continuous(expand = expansion(mult = 0.1)) +
    scale_y_continuous(expand = expansion(mult = 0.1)) +
    coord_cartesian(clip = "off") +
    labs(
      x = NULL, y = NULL, title = title,
      subtitle = sprintf(
        "%d proteins passing %d or more routes; edges: STRING score ≥ %d",
        nrow(shown), limits$community_routes, limits$string_score
      ),
      caption = if (length(alone)) {
        sprintf("%d proteins without a STRING partner (listed in the workbook)", length(alone))
      }
    ) +
    theme_figure() +
    theme(
      axis.text = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), legend.box = "vertical",
      plot.caption = element_text(size = 5.2, hjust = 0, colour = "grey30")
    )
  list(
    plot = plot,
    data = shown |>
      mutate(has_partner = gene %in% linked$gene) |>
      select(gene, logFC, P.Value, lfsr, routes, statistics:pathway, community, name, has_partner)
  )
}

# Every mouse for one sex's candidates passing CANDIDATE_DOT_ROUTES or more routes, most routes
# first and smallest p within: each protein's log2 abundance minus its own sex's vehicle mean, all
# four groups, open circles where the protein was not detected in that mouse. cands comes from
# 06_Candidates, fit from 01_Fit and quant from 02_Quantification. Returns the plot and the
# candidates drawn.
candidate_dots <- function(cands, fit, quant, sex, title, ncol = 7) {
  shown <- cands |>
    filter(sex == !!sex, routes >= CANDIDATE_DOT_ROUTES) |>
    arrange(desc(routes), P.Value)
  targets <- mutate(as_tibble(fit$targets), group = factor(group, FIGURE_GROUPS))
  values <- fit$E_norm$full[shown$feature, targets$sample_id, drop = FALSE]
  observations <- quant$proteins$other$n.observations
  detected <- observations[shown$feature, targets$sample_id, drop = FALSE] > 0
  genes <- set_names(
    sprintf("%s (%d)", protein_label(shown$gene, shown$feature), shown$routes), shown$feature
  )
  points <- as_tibble(centre_on_vehicle(values, targets), rownames = "feature") |>
    pivot_longer(-feature, names_to = "sample_id", values_to = "value") |>
    left_join(select(targets, sample_id, group), by = "sample_id") |>
    mutate(
      seen = detected[cbind(feature, sample_id)],
      protein = factor(genes[feature], unique(genes))
    )
  plot <- ggplot(points, aes(group, value, colour = group)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
    stat_summary(
      fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.25, colour = "grey30"
    ) +
    geom_point(aes(shape = seen),
      position = position_jitter(width = 0.12, seed = SEED), size = 0.9
    ) +
    facet_wrap(~protein, ncol = ncol, scales = "free_y") +
    scale_colour_manual(values = GROUP_COLOURS, labels = GROUP_LABELS, name = NULL) +
    scale_shape_manual(
      values = c(`TRUE` = 16, `FALSE` = 1),
      labels = c(`TRUE` = "detected", `FALSE` = "not detected"),
      name = NULL
    ) +
    scale_x_discrete(labels = \(g) str_wrap(GROUP_SHORT[g], 3)) +
    labs(
      x = NULL, y = "log2 vs own-sex vehicle", title = title,
      subtitle = sprintf(
        "Proteins passing %d or more routes (routes in brackets); bar: group mean",
        CANDIDATE_DOT_ROUTES
      )
    ) +
    theme_figure() +
    theme(
      axis.text.x = element_text(size = 4.5, lineheight = 0.9),
      strip.text = element_text(size = 5.5), panel.spacing = unit(1.2, "mm")
    ) +
    guides(colour = guide_legend(nrow = 1, order = 1), shape = guide_legend(order = 2))
  list(
    plot = plot,
    data = select(
      shown, gene, logFC, P.Value, adj.P.Val, lfsr, auc, routes, detected_veh, detected_pas
    )
  )
}

# Which of the six routes each candidate passed: one row per protein passing community_routes or
# more routes, most routes first and smallest p within, a filled tile per route passed, the log2
# FC as a square point (so the route tiles keep the fill scale to themselves) and the STRING
# neighbourhood's name to the right. Returns the plot and the candidates drawn.
candidate_routes <- function(cand, sex, title) {
  route_labels <- c(
    statistics = paste("p ≤", cand$thresholds$p), every_mouse = "Every mouse",
    both_sexes = "Both sexes", phenotype = "Phenotype", hub = "Module hub", pathway = "Leading edge"
  )
  shown <- cand$candidates |>
    filter(sex == !!sex, routes >= cand$thresholds$community_routes) |>
    arrange(desc(routes), P.Value) |>
    mutate(
      row = protein_label(gene, feature),
      row = factor(row, rev(row)),
      group = if_else(is.na(name), "", set_label(name, 46))
    )
  fc_column <- length(route_labels) + 1
  tiles <- shown |>
    select(row, all_of(names(route_labels))) |>
    pivot_longer(-row, names_to = "route", values_to = "passed") |>
    mutate(route = factor(route_labels[route], c(route_labels, "log2 FC")))
  limit <- max(abs(shown$logFC))
  plot <- ggplot(tiles, aes(route, row)) +
    geom_tile(aes(fill = passed), colour = "white", linewidth = 0.4) +
    scale_fill_manual(
      values = c(`TRUE` = "grey25", `FALSE` = "grey92"), guide = "none"
    ) +
    geom_point(aes(x = fc_column, y = row, colour = logFC), shown, shape = 15, size = 2.2) +
    scale_colour_gradient2(
      low = DIRECTION_COLOURS[["Down"]], high = DIRECTION_COLOURS[["Up"]],
      limits = c(-limit, limit), name = "log2 FC", guide = colour_bar()
    ) +
    geom_text(aes(x = fc_column + 0.7, label = group), shown,
      hjust = 0, size = 1.6, colour = "grey30"
    ) +
    scale_x_discrete(position = "top", drop = FALSE, expand = expansion(add = c(0.6, 9.5))) +
    coord_cartesian(clip = "off") +
    labs(
      x = NULL, y = NULL, title = title,
      subtitle = "Filled: route passed; right: STRING neighbourhood"
    ) +
    theme_figure() +
    theme(
      panel.grid.major = element_blank(),
      axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 5.5),
      axis.text.y = element_text(size = 5.2, face = "italic")
    )
  list(
    plot = plot,
    data = select(shown, gene, routes, all_of(names(route_labels)), logFC, P.Value, lfsr, name)
  )
}
