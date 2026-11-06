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
MODEL_LABELS <- c(
  pooled = "Pooled", pooled_mito = "Pooled + pellet MitoCarta share",
  treatment = "Within treatment"
)
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
# Dot area stops growing at this -log10 FDR or p.
LOG_CAP <- 4

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

# Gene symbols, with the accession added where two protein groups share one; a protein may appear
# on several rows.
protein_label <- function(gene, feature) {
  groups <- distinct(tibble(gene, feature))
  shared <- gene %in% groups$gene[duplicated(groups$gene)]
  if_else(shared, paste(gene, feature), gene)
}

# p as the labels print it, without scientific notation.
format_p <- \(p) if_else(p < 0.001, "p < 0.001", sprintf("p = %.2g", p))

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

# Mice as points, fill by group and shape by treatment (circle vehicle, triangle PAS), with one
# key that shows both; the data need group (as FIGURE_GROUPS codes). position goes to geom_point().
group_points <- function(position = "identity") {
  list(
    geom_point(aes(fill = group, shape = str_sub(group, 1, 3)),
      position = position, size = 1.5, colour = "grey15", stroke = 0.25
    ),
    scale_fill_manual(
      values = GROUP_COLOURS, labels = GROUP_LABELS, name = NULL,
      guide = guide_legend(
        nrow = 2, byrow = TRUE, order = 1, override.aes = list(shape = c(21, 24, 21, 24))
      )
    ),
    scale_shape_manual(values = c(VEH = 21, PAS = 24), guide = "none")
  )
}

# A main figure's folder part (a_script, b_reports); a supplementary figure's files, named S and
# its number, sit in that part's supp/.
figure_dir <- function(figure, part, name) {
  dir <- here("04_Figures", figure, part)
  if (str_detect(name, "^S[0-9]")) file.path(dir, "supp") else dir
}

# A panel or composite at its printed size into b_reports/ (supp/ for a supplement), as PDF
# (cairo, which carries Π and ≤) and as a 600 dpi PNG; panels are saved at the size they take in
# the composite, so it stitches panels that already look final.
save_figure <- function(plot, figure, name, width_mm, height_mm) {
  dir <- figure_dir(figure, "b_reports", name)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  save <- \(ext, ...) {
    ggsave(file.path(dir, paste0(name, ext)), plot,
      width = width_mm, height = height_mm, units = "mm", ...
    )
  }
  save(".pdf", device = grDevices::cairo_pdf)
  save(".png", device = ragg::agg_png, dpi = 600)
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
# contrast_strips() key. Panels have a fixed size, 5 mm a column and 3 mm a row, a block at least
# min_block_mm tall so a short block's strip label fits. Dot area stops growing at LOG_CAP.
pathway_grid <- function(cells, rows, columns, selected_label, title = NULL, subtitle = NULL,
                         tones = NULL, fill_name = "NES", fill_limit = NULL,
                         fill_breaks = c(-2, 0, 2), size_name = "FDR", outline_shape = 22,
                         min_block_mm = 0) {
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
        c(if (length(lead)) 3.5, 5 * table(columns$column_block)[column_blocks]), "mm"
      ),
      rows = unit(block_mm(rows, min_block_mm), "mm")
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
    scale_size(
      range = c(0.6, 2.8), limits = c(0, LOG_CAP), breaks = -log10(c(0.1, 0.05, 0.01, 10^-LOG_CAP)),
      labels = c("0.1", "0.05", "0.01", paste("≤", format(10^-LOG_CAP, scientific = FALSE))),
      name = size_name
    ) +
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
      axis.text.x = element_text(size = 5.5, angle = 90, hjust = 1, vjust = 0.5),
      strip.placement = "outside", strip.text = element_text(size = 5.5),
      strip.text.y.left = element_text(angle = 90), panel.spacing = unit(1, "mm"),
      legend.position = "right", legend.box.just = "left", legend.key.size = unit(2.6, "mm"),
      legend.text = element_text(size = 5.5), legend.title = element_text(size = 6)
    )
}

# A pathway_grid()'s blocks in mm, top to bottom; their sum and the 1 mm between blocks are the
# grid's body, to which a figure adds its title, column labels and keys.
block_mm <- function(rows, min_block_mm = 0) {
  blocks <- levels(droplevels(rows$block))
  pmax(3 * table(rows$block)[blocks], min_block_mm)
}
grid_body_mm <- \(rows, min_block_mm = 0) {
  sum(block_mm(rows, min_block_mm)) + nlevels(droplevels(rows$block)) - 1
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
      column_block = names(SEX_TRAIT_TONES)[match(sex, SEX_NAMES)],
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
count_bars <- function(counts, pairs, rows, width_mm = 18) {
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
    labs(
      x = sprintf("Distinct significant sets\n(FDR ≤ %s, after overlap grouping)", ALPHA),
      y = NULL
    ) +
    theme_figure() +
    theme(
      axis.text.y = element_blank(), panel.grid.major.y = element_blank(),
      strip.text = element_blank(), panel.spacing = unit(1, "mm"),
      axis.text.x = element_text(size = 5), axis.title.x = element_text(size = 5.5),
      plot.margin = margin(l = 22, r = 4)
    )
}

# Display names for pathway sets: the collection prefix dropped, title case, shortened to width
# characters when a width is given.
set_label <- function(term, width = NULL) {
  label <- enrichVolcano::clean_label(term, width = 200)
  if (is.null(width)) label else str_trunc(label, width)
}

# A figure's panel scripts, each run in its own environment, as YvO's stitchers do, and each plot
# given its letter when there is more than one. A panel script returns list(plot, data).
run_panels <- function(figure, panels) {
  map(set_names(panels), \(name) {
    path <- file.path(figure_dir(figure, "a_script", name), paste0(name, ".R"))
    panel <- source(path, local = new.env())$value
    if (length(panels) > 1) panel$plot <- panel$plot + labs(tag = letters[match(name, panels)])
    panel
  })
}

# A figure's workbook, c_data/<name>.xlsx: the table each panel draws, one sheet per panel (a
# panel whose data is a list of tables gives one sheet per table, named <panel>_<table>).
write_figure_workbook <- function(panels, figure, name) {
  sheets <- list_flatten(map(panels, "data"), name_spec = "{outer}_{inner}")
  file <- here("04_Figures", figure, "c_data", paste0(name, ".xlsx"))
  dir.create(dirname(file), showWarnings = FALSE)
  writexl::write_xlsx(sheets, file)
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
      section = names(SEX_TRAIT_TONES)[match(sex, SEX_NAMES)],
      column = TRAIT_SHORT[phenotype], r_equiv, pearson_r = r, P.Value, adj.P.Val
    )
  bind_rows(contrast_cells, trait_cells) |>
    mutate(
      section = factor(section, c("Contrast", names(SEX_TRAIT_TONES))),
      column = factor(column, c(CONTRAST_SHORT, TRAIT_SHORT))
    )
}
