# The steps every arc runs on its own slice of the trunk. Every model, test and score is computed
# once in 00_Trunk (07_Fit, 08_Enrichment, 09_Classification); an arc step reads those results,
# keeps its contrast or sex, and writes one workbook (every feature, every statistic), one report
# PDF (overview pages, then a grid of every feature at p <= 0.05) and the rds its figures read.
# Sourced by the arc scripts; it sources R/helpers.R.

source(here::here("R", "helpers.R"))

GRID_PER_PAGE <- 20

trunk_path <- \(step, file) here("00_Trunk", step, "c_data", file)

arc_paths <- \(arc, step) stage_paths(ARCS[[arc]]$dir, step)

# The mice an arc looks at: one sex for the PAS arcs, all twenty for the sex difference.
arc_mice <- function(arc, targets) {
  sex <- ARCS[[arc]]$sex
  if (is.na(sex)) targets$sample_id else targets$sample_id[targets$sex == sex]
}

arc_colour <- \(arc) unname(CONTRAST_COLOURS[CONTRAST_LABELS[ARCS[[arc]]$contrast]])

# The strongest tier a protein reaches: FDR, then Π, then nominal p.
tier_of <- function(p, fdr, pi) {
  factor(
    case_when(
      fdr <= ALPHA ~ "FDR ≤ 0.05", pi <= ALPHA ~ "Π ≤ 0.05",
      p <= ALPHA ~ "p ≤ 0.05",
      .default = "p > 0.05"
    ),
    c("p > 0.05", "p ≤ 0.05", "Π ≤ 0.05", "FDR ≤ 0.05")
  )
}
tier_colours <- function(arc) {
  colour <- arc_colour(arc)
  c(
    "p > 0.05" = "grey82", "p ≤ 0.05" = grDevices::colorRampPalette(c(colour, "white"))(101)[56],
    "Π ≤ 0.05" = colour, "FDR ≤ 0.05" = "black"
  )
}

# Splits long data into report pages of `per_page` facets, keeping the order of `strip`; a short
# last page is padded with blank facets so every page keeps the same grid.
paginate <- function(data, draw, per_page = GRID_PER_PAGE) {
  strips <- levels(data$strip)
  if (length(strips) == 0) {
    return(list())
  }
  pad <- strrep(" ", seq_len(per_page))
  split(strips, ceiling(seq_along(strips) / per_page)) |>
    map(\(s) {
      levels <- c(s, pad[seq_len(per_page - length(s))])
      draw(filter(data, strip %in% s) |> mutate(strip = factor(strip, levels)))
    })
}

# Facet labels: a feature's name and its statistics, in the order given.
strip_factor <- function(label) factor(label, unique(label))
short <- \(x, width = 34) str_trunc(enrichVolcano::clean_label(x, width = 200), width)

# Per-mouse values of many features, one facet each, grouped and coloured by group.
mouse_grid <- function(long, title, subtitle, y = "log2 abundance") {
  paginate(long, \(d) {
    ggplot(d, aes(group, value, colour = group, shape = sex)) +
      stat_summary(
        fun = median, geom = "crossbar", width = 0.5, linewidth = 0.25, colour = "grey40"
      ) +
      geom_point(position = position_jitter(width = 0.12, seed = SEED), size = 1.1) +
      facet_wrap(~strip, scales = "free_y", ncol = 5, nrow = 4, drop = FALSE) +
      scale_colour_manual(values = GROUP_COLOURS, guide = "none") +
      scale_shape_manual(values = SEX_SHAPES, guide = "none") +
      labs(x = NULL, y = y, title = title, subtitle = subtitle) +
      theme_figure() +
      theme(strip.text = element_text(size = 5, face = "plain"), axis.text = element_text(size = 5))
  })
}
long_values <- function(m, features, strips, targets) {
  as_tibble(m[features, targets$sample_id, drop = FALSE], rownames = "feature") |>
    mutate(strip = strips) |>
    pivot_longer(-c(feature, strip), names_to = "sample_id") |>
    left_join(select(targets, sample_id, group, sex, treatment), by = "sample_id") |>
    mutate(group = factor(group, FIGURE_GROUPS), strip = factor(strip, unique(strips)))
}

# Proteins: the arc's contrast sliced from the shared fit, its tiers and Π lists, and the STRING
# connectivity of those lists from 08_Enrichment (descriptive).
arc_proteins <- function(arc) {
  out <- arc_paths(arc, "01_Proteins")
  inputs <- c(
    fit = trunk_path("07_Fit", "fit.rds"),
    enrichment = trunk_path("08_Enrichment", "enrichment.rds")
  )
  fit <- readRDS(inputs[["fit"]])
  contrast <- ARCS[[arc]]$contrast
  label <- ARCS[[arc]]$label
  targets <- filter(fit$targets, sample_id %in% arc_mice(arc, fit$targets))

  results <- fit$protein_contrasts |>
    filter(contrast == !!contrast) |>
    mutate(
      tier = tier_of(P.Value, adj.P.Val, pi_score), direction = if_else(logFC > 0, "up", "down")
    ) |>
    arrange(universe, P.Value)
  summary <- results |>
    reframe(tier_counts(P.Value, adj.P.Val, pi_score), .by = c(universe, direction))
  pi_lists <- results |>
    filter(universe == "full", pi_score <= ALPHA, !is.na(gene), gene != "") |>
    select(direction, feature, gene, logFC, P.Value, adj.P.Val, pi_score, r_group) |>
    arrange(direction, pi_score)
  string <- filter(readRDS(inputs[["enrichment"]])$string, contrast == !!contrast)

  volcano <- results |>
    mutate(universe = if_else(universe == "full", "Full proteome", "Mitochondrial universe")) |>
    ggplot(aes(logFC, -log10(P.Value), colour = tier)) +
    geom_point(size = 0.6) +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", linewidth = 0.3) +
    ggrepel::geom_text_repel(
      data = \(d) filter(d, adj.P.Val <= ALPHA | rank(P.Value) <= 10, .by = universe),
      aes(label = gene), size = 2, colour = "black", max.overlaps = Inf,
      min.segment.length = 0, seed = SEED
    ) +
    facet_wrap(~universe) +
    scale_colour_manual(values = tier_colours(arc), name = NULL, drop = FALSE) +
    labs(
      x = "log2 fold change", y = "-log10 p",
      title = str_glue("{label}: {CONTRAST_LABELS[contrast]}, every protein"),
      subtitle = "limpa dpcDE, BH within universe; labelled: FDR ≤ 0.05 and the 10 smallest p"
    ) +
    theme_figure() +
    theme(legend.position = "bottom")
  tiers <- c("p_le_0.05", "pi_le_0.05", "fdr_le_0.05")
  counts <- summary |>
    pivot_longer(all_of(tiers), names_to = "tier", values_to = "n") |>
    mutate(
      tier = factor(tier, tiers, c("p ≤ 0.05", "Π ≤ 0.05", "FDR ≤ 0.05")),
      universe = if_else(universe == "full", "Full proteome", "Mitochondrial universe"),
      n = if_else(direction == "down", -n, n)
    ) |>
    ggplot(aes(n, tier, fill = direction)) +
    geom_col(width = 0.7) +
    geom_vline(
      aes(xintercept = -expected_by_chance),
      data = \(d) filter(d, tier == "p ≤ 0.05", direction == "down"),
      linetype = "dashed", linewidth = 0.3
    ) +
    geom_vline(
      aes(xintercept = expected_by_chance),
      data = \(d) filter(d, tier == "p ≤ 0.05", direction == "up"),
      linetype = "dashed", linewidth = 0.3
    ) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_text(aes(label = abs(n), hjust = if_else(n < 0, 1.2, -0.2)), size = 2.2) +
    facet_wrap(~universe, scales = "free_x") +
    scale_fill_manual(values = c(down = "#2166AC", up = "#B2182B"), name = NULL) +
    scale_x_continuous(labels = abs, expand = expansion(mult = 0.2)) +
    labs(
      x = "Proteins (down | up)", y = NULL, title = "Proteins per tier",
      subtitle = paste(
        "Dashed: proteins at p ≤ 0.05 expected by chance in each direction;",
        "Π can include p > 0.05"
      )
    ) +
    theme_figure()
  hits <- filter(results, universe == "full", P.Value <= ALPHA)
  heatmap <- group_heatmap(fit$E_norm$full[hits$feature, targets$sample_id, drop = FALSE], targets,
    title = str_glue("All {nrow(hits)} proteins at p ≤ 0.05, per mouse"),
    subtitle = "Rows z-scored across these mice, clustered on 1 - Pearson r (display only)"
  )
  interaction <- is.na(ARCS[[arc]]$sex)
  strips <- with(hits, sprintf(
    "%s\nlog2FC %.2f%s\np %.2g, FDR %.2g, Π %.2g", gene, logFC,
    if (interaction) "" else sprintf(", r %.2f", r_group), P.Value, adj.P.Val, pi_score
  ))
  grids <- mouse_grid(long_values(fit$E_norm$full, hits$feature, strips, targets),
    title = str_glue("{label}: every protein at p ≤ 0.05, in order of p"),
    subtitle = if (interaction) {
      paste(
        "log2 abundance per mouse, bar at the group median;",
        "log2FC: female PAS response minus male PAS response"
      )
    } else {
      "log2 abundance per mouse, bar at the group median; r: Pearson r with PAS (PAS = 1)"
    }
  )
  plot_index <- write_report(file.path(out$reports, "01_proteins.pdf"), c(
    list(
      "volcanoes, every protein" = volcano, "proteins per tier" = counts,
      "heatmap of proteins at p ≤ 0.05" = heatmap
    ),
    set_names(grids, sprintf("per-mouse grid, page %d", seq_along(grids)))
  ))

  saveRDS(
    list(results = results, summary = summary, pi_lists = pi_lists, string = string),
    file.path(out$data, "proteins.rds")
  )
  write_workbook(
    list(
      summary = summary, results = select(results, -genes), pi_lists = pi_lists, string = string,
      plot_index = plot_index
    ),
    c(
      "Per universe and direction: p ≤ 0.05, chance, Π ≤ 0.05, FDR tiers.",
      paste(
        "Every protein in both universes, by p:",
        "log2 FC with 95% CI, t, p, BH FDR, Π, r with PAS, tier."
      ),
      "Full-universe proteins at Π ≤ 0.05, by direction.",
      "STRING v12 connectivity of each list against the detected proteome (descriptive).",
      "Which PDF page shows what."
    ),
    out$workbook, inputs
  )
  summary
}

# EnrichmentMap similarity (mean of Jaccard and overlap coefficient; Merico et al. 2010), the
# measure enrichVolcano's collapse uses, so the tree and the collapse agree.
enrichmentmap_similarity <- function(sets) {
  outer(names(sets), names(sets), Vectorize(\(a, b) {
    shared <- length(intersect(sets[[a]], sets[[b]]))
    mean(c(
      shared / length(union(sets[[a]], sets[[b]])),
      shared / min(length(sets[[a]]), length(sets[[b]]))
    ))
  })) |>
    `dimnames<-`(list(names(sets), names(sets)))
}

# fgsea running-score plots, twelve per page.
enrichment_grid <- function(terms, sets, ranks, title) {
  if (nrow(terms) == 0) {
    return(list())
  }
  plots <- pmap(terms, \(term, NES, p, padj, size, universe, score_p, ...) {
    fgsea::plotEnrichment(sets[[universe]][[term]], ranks[[universe]], ticksSize = 0.1) +
      labs(
        x = NULL, y = NULL,
        title = sprintf(
          "%s\nNES %.2f, p %.2g, FDR %.2g, n %d\nsingscore p %.2g", short(term, 40), NES, p, padj,
          size, score_p
        )
      ) +
      theme_figure() +
      theme(
        plot.title = element_text(size = 5, face = "plain"), axis.text = element_text(size = 4.5)
      )
  })
  split(plots, ceiling(seq_along(plots) / 12)) |>
    map(\(page) {
      patchwork::wrap_plots(
        c(page, rep(list(patchwork::plot_spacer()), 12 - length(page))),
        ncol = 4
      ) +
        patchwork::plot_annotation(
          title = title, theme = theme(plot.title = element_text(size = 7.5, face = "bold"))
        )
    })
}

# Pathways: the arc's contrast sliced from 08_Enrichment. fgsea on GO:BP and Reactome in the full
# proteome and MitoCarta in the mitochondrial universe, BH within collection over every tested set,
# each beside the same set's singscore contrast; GO Slim compartments and the literature list
# alongside. GO Slim generic BP families (every slim ancestor a term reaches) organise the GO:BP
# terms; they never choose what is tested. Rings take the 12 strongest terms per collection that
# the EnrichmentMap flag keeps; the tree takes the 20 strongest full-proteome terms. Display only.
arc_pathways <- function(arc) {
  out <- arc_paths(arc, "02_Pathways")
  inputs <- c(
    enrichment = trunk_path("08_Enrichment", "enrichment.rds"),
    gene_sets = trunk_path("04_Gene_Sets", "gene_sets.rds")
  )
  enr <- readRDS(inputs[["enrichment"]])
  gs <- readRDS(inputs[["gene_sets"]])
  contrast <- ARCS[[arc]]$contrast
  label <- ARCS[[arc]]$label
  COLLECTIONS <- c(`GO:BP` = "full", Reactome = "full", MitoCarta = "mito")

  results <- filter(enr$results, contrast == !!contrast) |>
    select(-contrast) |>
    arrange(factor(collection, names(COLLECTIONS)), p)
  summary <- filter(enr$summary, contrast == !!contrast) |>
    select(-contrast)
  compartments <- filter(enr$compartments, contrast == !!contrast) |>
    select(-contrast)
  literature <- filter(enr$literature, contrast == !!contrast) |>
    select(-contrast)
  families <- gs$goslim |>
    transmute(term, family = str_split(slim, "; ")) |>
    unnest(family) |>
    filter(family != "")
  family_summary <- results |>
    filter(collection == "GO:BP") |>
    inner_join(families, by = "term", relationship = "many-to-many") |>
    summarise(
      terms = n(), up_p = sum(p <= ALPHA & NES > 0), down_p = sum(p <= ALPHA & NES < 0),
      up_fdr = sum(padj <= ALPHA & NES > 0), down_fdr = sum(padj <= ALPHA & NES < 0),
      .by = family
    ) |>
    arrange(desc(up_p + down_p))

  display <- filter(results, !(collection == "GO:BP" & term %in% gs$slim_terms))
  ring_terms <- display |>
    filter(dedup_status == "kept") |>
    slice_min(p, n = 12, by = collection, with_ties = FALSE)
  tree_terms <- display |>
    filter(universe == "full") |>
    slice_min(p, n = 20, with_ties = FALSE)
  tree <- hclust(as.dist(1 - enrichmentmap_similarity(enr$sets$full[tree_terms$term])), "average")
  tree_terms <- tree_terms |>
    mutate(cluster = cutree(tree, h = 1 - 0.375)[term]) |>
    mutate(cluster_name = term[which.min(p)], .by = cluster)

  rings <- imap(COLLECTIONS, \(u, db) {
    terms <- filter(ring_terms, collection == db)
    da <- enr$da[[u]]
    enrichVolcano::plot_volcano_ring(da[da$contrast == contrast, ], enr$objects[[contrast]][[u]],
      contrast = contrast, databases = db, terms = terms$term,
      labels = set_names(
        paste0(
          enrichVolcano::clean_label(terms$term, width = 22),
          if_else(terms$padj <= ALPHA, " *", "")
        ),
        terms$term
      ),
      title = str_glue("{label}: {db}"),
      subtitle = paste0(
        "12 strongest terms the EnrichmentMap flag keeps; * FDR ≤ 0.05",
        if (is.na(ARCS[[arc]]$sex)) "\nup: PAS raised the set more in females than in males" else ""
      ),
      p_threshold = ALPHA, show_counts = FALSE, score_limits = c(-3, 3), arc_span = "auto",
      theme = enrichVolcano::plot_theme(base_size = 7, base_family = "Helvetica")
    )
  })
  family_page <- family_summary |>
    filter(up_p + down_p > 0) |>
    pivot_longer(c(up_p, down_p), names_to = "direction", values_to = "n") |>
    mutate(
      n = if_else(direction == "down_p", -n, n),
      family = reorder(str_to_sentence(family), abs(n), sum)
    ) |>
    ggplot(aes(n, family, fill = direction)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_fill_manual(
      values = c(down_p = "#2166AC", up_p = "#B2182B"), labels = c("down", "up"), name = NULL
    ) +
    scale_x_continuous(labels = abs) +
    labs(
      x = "GO:BP terms at p ≤ 0.05 (down | up)", y = NULL,
      title = str_glue("{label}: signal by GO Slim family"),
      subtitle = paste(
        "fgsea; a term counts in every family it belongs to;",
        "every family with a term at p ≤ 0.05"
      )
    ) +
    theme_figure()
  tree_page <- tree_terms |>
    mutate(term = factor(term, tree$labels[tree$order])) |>
    ggplot(aes(NES, term)) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_point(aes(size = -log10(p), fill = NES), shape = 21, colour = "grey60", stroke = 0.2) +
    geom_point(
      data = \(d) filter(d, padj <= ALPHA), aes(size = -log10(p)), shape = 21, fill = NA,
      colour = "black", stroke = 0.7
    ) +
    scale_y_discrete(labels = \(x) enrichVolcano::clean_label(x, width = 60)) +
    scale_fill_gradient2(
      low = "#2166AC", high = "#B2182B", limits = c(-3, 3), oob = scales::squish
    ) +
    labs(
      x = "NES", y = NULL, title = "20 strongest GO:BP and Reactome terms, ordered by shared genes",
      subtitle = "EnrichmentMap similarity, average linkage; black ring FDR ≤ 0.05"
    ) +
    theme_figure()
  compartment_page <- compartments |>
    mutate(compartment = reorder(enrichVolcano::clean_label(set, width = 40), NES)) |>
    ggplot(aes(NES, compartment, fill = padj <= ALPHA)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_text(aes(label = sprintf("p %.2g", p), hjust = if_else(NES < 0, 1.1, -0.1)), size = 1.8) +
    scale_fill_manual(
      values = c(`TRUE` = arc_colour(arc), `FALSE` = "grey70"),
      labels = c("FDR > 0.05", "FDR ≤ 0.05"), name = NULL
    ) +
    scale_x_continuous(expand = expansion(mult = 0.25)) +
    labs(
      x = "NES", y = NULL, title = "GO Slim cellular compartments",
      subtitle = "fgsea, full proteome"
    ) +
    theme_figure()
  grids <- enrichment_grid(filter(results, p <= ALPHA), enr$sets, enr$ranks[[contrast]],
    title = str_glue("{label}: every term at p ≤ 0.05, in order of p (fgsea running score)")
  )
  plot_index <- write_report(file.path(out$reports, "02_pathways.pdf"), c(
    set_names(rings, paste("ring,", names(rings))),
    list(
      "signal by GO Slim family" = family_page,
      "strongest full-proteome terms by shared genes" = tree_page, "compartments" = compartment_page
    ),
    set_names(grids, sprintf("running-score grid, page %d", seq_along(grids)))
  ))

  saveRDS(
    list(
      compartments = compartments, results = results, ring_terms = ring_terms, tree = tree,
      tree_terms = tree_terms, families = family_summary, literature = literature
    ),
    file.path(out$data, "pathways.rds")
  )
  write_workbook(
    list(
      summary = summary, fgsea = results, families = family_summary, compartments = compartments,
      ring_terms = select(ring_terms, universe, collection, term, NES, p, padj),
      tree_terms = select(tree_terms, term, NES, p, padj, cluster, cluster_name),
      literature = literature, plot_index = plot_index
    ),
    c(
      "Per universe and collection: p ≤ 0.05, chance, FDR tiers.",
      paste(
        "Every tested set, by p: NES, p, BH FDR within collection, size, redundancy flag,",
        "singscore",
        "contrast (log2 FC, p, FDR, r with group), GO Slim families, leading edge."
      ),
      "Per GO Slim family: GO:BP terms tested and at p and FDR ≤ 0.05 by direction.",
      "GO Slim cellular compartments: fgsea, full proteome.",
      "Ring terms: 12 strongest per collection that the EnrichmentMap flag keeps.",
      paste(
        "Tree terms: 20 strongest GO:BP and Reactome terms,",
        "their cluster at EnrichmentMap 0.375 and its lead term."
      ),
      "Literature-informed and ROS sets, fgsea in the full proteome, BH across the list.",
      "Which PDF page shows what."
    ),
    out$workbook, inputs
  )
  summary
}

# Scatter of feature against phenotype per facet: a line within vehicle, a line within PAS and a
# dashed line through all the sex's mice, the association that is tested.
link_grid <- function(long, title, subtitle, ncol = 5, nrow = 4) {
  paginate(long, per_page = ncol * nrow, draw = \(d) {
    ggplot(d, aes(phenotype_value, value)) +
      geom_smooth(
        aes(group = 1),
        method = "lm", formula = y ~ x, se = FALSE, colour = "grey20", linetype = "dashed",
        linewidth = 0.4
      ) +
      geom_smooth(
        aes(colour = group),
        method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.4
      ) +
      geom_point(aes(colour = group), size = 1.1) +
      facet_wrap(~strip, scales = "free", ncol = ncol, nrow = nrow, drop = FALSE) +
      scale_colour_manual(values = GROUP_COLOURS, name = NULL) +
      labs(x = "Phenotype", y = "Feature value", title = title, subtitle = subtitle) +
      theme_figure() +
      theme(
        strip.text = element_text(size = 4.6, face = "plain"), axis.text = element_text(size = 4.5),
        legend.position = "bottom"
      )
  })
}

# The protein with the smallest association p whose slope keeps its sign and p <= 0.05 when any
# one mouse is left out (ordinary least squares within the sex), so an example is not one mouse's
# doing.
robust_example <- function(links, E, x, mice, max_tries = 50) {
  candidates <- slice_min(links, P.Value, n = max_tries, with_ties = FALSE)
  for (i in seq_len(nrow(candidates))) {
    y <- E[candidates$feature[i], mice]
    holds <- map_lgl(seq_along(mice), \(k) {
      f <- summary(lm(y[-k] ~ x[-k]))$coefficients
      sign(f[2, 1]) == sign(candidates$estimate[i]) && f[2, 4] <= ALPHA
    })
    if (all(holds)) {
      return(candidates[i, ])
    }
  }
  candidates[0, ]
}

LEVEL_LABELS <- c(
  protein = "Proteins", pathway = "Pathways", module = "Modules", complex = "Complexes",
  cluster = "Clusters"
)

# Associations, two different questions kept apart. Phenotype effects: what PAS did to each
# phenotype in this sex (% change against vehicle, a treatment effect). Phenotype associations:
# which proteins, pathway scores (GO:BP, Reactome, MitoCarta), modules, complexes and clusters move
# with each phenotype across this sex's ten PAS and vehicle mice (07_Fit; r pooled and within each
# group beside every test). An association is co-movement, not a PAS effect. Pathway patterns also
# come from fgsea on the protein association t (08_Enrichment); that test permutes genes, so its FDR
# runs optimistic against the mouse-level singscore association of the same set.
arc_associations <- function(arc) {
  out <- arc_paths(arc, "03_Associations")
  inputs <- c(
    fit = trunk_path("07_Fit", "fit.rds"),
    enrichment = trunk_path("08_Enrichment", "enrichment.rds")
  )
  fit <- readRDS(inputs[["fit"]])
  sex <- ARCS[[arc]]$sex
  label <- ARCS[[arc]]$label
  targets <- filter(fit$targets, sample_id %in% arc_mice(arc, fit$targets))
  mice <- targets$sample_id
  phenotype_of <- \(p) fit$phenotypes[[p]][match(mice, fit$targets$sample_id)]

  effects <- filter(fit$phenotype_effects, sex == !!sex)
  links <- bind_rows(mutate(fit$protein_links, level = "protein"), fit$score_links) |>
    filter(sex == arc) |>
    mutate(label = coalesce(gene, feature), .after = feature)
  link_summary <- links |>
    reframe(tier_counts(P.Value, adj.P.Val), .by = c(level, universe, phenotype, adjusted)) |>
    select(-pi_le_0.05)
  links <- filter(links, !adjusted) |>
    arrange(level, phenotype, P.Value)
  r_group <- bind_rows(
    mutate(fit$protein_contrasts, level = "protein"), fit$score_contrasts
  ) |>
    filter(contrast == ARCS[[arc]]$contrast) |>
    select(level, universe, feature, r_group, group_p = P.Value, group_fdr = adj.P.Val)
  fgsea <- readRDS(inputs[["enrichment"]])$association_fgsea |>
    filter(sex == arc) |>
    select(-sex)

  examples <- map(PHENOTYPES, \(p) {
    links |>
      filter(level == "protein", universe == "full", phenotype == p) |>
      robust_example(fit$E_norm$full, phenotype_of(p), mice)
  }) |>
    list_rbind()

  effects_page <- phenotype_effect_plot(effects,
    title = str_glue("{label}: phenotype effects of PAS"),
    subtitle = "% change against vehicle, Welch t on log values, n = 5 v 5"
  )
  grid_page <- link_summary |>
    filter(!adjusted) |>
    mutate(
      level = factor(
        paste0(LEVEL_LABELS[level], if_else(universe == "mito" & level == "protein", ", mito", "")),
        c("Proteins", "Proteins, mito", LEVEL_LABELS[-1])
      ),
      trait = factor(PHENOTYPE_LABELS[phenotype], rev(PHENOTYPE_LABELS))
    ) |>
    summarise(
      across(c(tests, p_le_0.05, expected_by_chance, fdr_le_0.05), sum),
      .by = c(level, trait)
    ) |>
    mutate(ratio = p_le_0.05 / expected_by_chance) |>
    ggplot(aes(level, trait, fill = ratio)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = str_glue("{p_le_0.05} / {tests}\nFDR {fdr_le_0.05}")), size = 1.9) +
    scale_fill_gradient(low = "white", high = "grey35", name = "p ≤ 0.05 /\nchance") +
    labs(
      x = NULL, y = NULL, title = str_glue("{label}: phenotype associations by level"),
      subtitle = paste(
        "Co-movement across PAS and vehicle mice of this sex (n = 10);",
        "p ≤ 0.05 of tests; BH within level and phenotype"
      )
    ) +
    theme_figure()
  check_page <- link_summary |>
    filter(level == "protein", universe == "full") |>
    select(phenotype, adjusted, p_le_0.05, fdr_le_0.05) |>
    pivot_longer(c(p_le_0.05, fdr_le_0.05), names_to = "count") |>
    mutate(
      trait = factor(PHENOTYPE_LABELS[phenotype], rev(PHENOTYPE_LABELS)),
      count = if_else(count == "p_le_0.05", "p ≤ 0.05", "FDR ≤ 0.05"),
      model = if_else(adjusted, "with pellet mito share", "without")
    ) |>
    ggplot(aes(value, trait, fill = model)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    facet_wrap(~count, scales = "free_x") +
    scale_fill_manual(
      values = c(without = "grey40", `with pellet mito share` = "grey75"), name = NULL
    ) +
    labs(
      x = "Proteins", y = NULL,
      title = "Check: protein associations with and without pellet mitochondrial share",
      subtitle = "Share centred within group, added to the association model"
    ) +
    theme_figure()
  group_rows <- links |>
    filter(level %in% c("module", "complex", "cluster"), universe == "full") |>
    select(level, feature, phenotype, r, adj.P.Val) |>
    left_join(r_group, by = c("level", "feature"), relationship = "many-to-one")
  heat <- bind_rows(
    distinct(group_rows, level, feature, r = r_group, fdr = group_fdr) |>
      mutate(column = "PAS vs vehicle"),
    transmute(group_rows, level, feature, r, fdr = adj.P.Val, column = PHENOTYPE_LABELS[phenotype])
  ) |>
    mutate(
      column = factor(column, c("PAS vs vehicle", PHENOTYPE_LABELS)),
      panel = if_else(column == "PAS vs vehicle", "Treatment effect", "Phenotype associations"),
      feature = str_trunc(enrichVolcano::clean_label(feature, width = 200), 45, side = "center"),
      level = LEVEL_LABELS[level]
    )
  heat_pages <- heat |>
    mutate(
      chunk = ceiling(match(feature, unique(feature)) / 32), .by = level
    ) |>
    mutate(level = factor(level, LEVEL_LABELS)) |>
    group_by(level, chunk)
  keys <- group_keys(heat_pages)
  r_heatmaps <- group_split(heat_pages) |>
    set_names(sprintf("%s, page %d", keys$level, keys$chunk)) |>
    map(\(d) {
      ggplot(d, aes(column, feature, fill = r)) +
        geom_tile(colour = "white") +
        geom_tile(
          data = \(x) filter(x, fdr <= ALPHA), fill = NA, colour = "black", linewidth = 0.5
        ) +
        geom_text(aes(label = sprintf("%.2f", r)), size = 1.7) +
        facet_grid(~panel, scales = "free_x", space = "free_x") +
        scale_fill_gradient2(low = "#2166AC", high = "#B2182B", limits = c(-1, 1)) +
        labs(
          x = NULL, y = NULL,
          title = str_glue("{label}: {d$level[1]}, treatment effect and phenotype associations"),
          subtitle = paste(
            "Pearson r across this sex's ten mice (treatment: PAS = 1);",
            "box: FDR ≤ 0.05 (limma, BH within level)"
          )
        ) +
        theme_figure() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
    })
  example_long <- examples |>
    rowwise() |>
    reframe(
      strip = sprintf(
        "%s · %s\nr %.2f (VEH %.2f, PAS %.2f)\np %.2g, FDR %.2g",
        PHENOTYPE_LABELS[phenotype], label, r,
        r_veh, r_pas, P.Value, adj.P.Val
      ),
      sample_id = mice, value = fit$E_norm$full[feature, mice],
      phenotype_value = phenotype_of(phenotype)
    ) |>
    left_join(select(targets, sample_id, group), by = "sample_id") |>
    mutate(group = factor(group, FIGURE_GROUPS), strip = strip_factor(strip))
  examples_page <- link_grid(example_long,
    title = str_glue("{label}: strongest robust protein per phenotype"),
    subtitle =
      "Smallest association p that keeps its sign and p ≤ 0.05 with any one mouse left out",
    ncol = 3, nrow = 2
  )

  score_matrix <- rbind(fit$scores$full, fit$scores$mito, fit$scores$module)
  feature_matrix <- \(level, universe) {
    if (level == "protein") fit$E_norm[[universe]] else score_matrix
  }
  hits <- links |>
    filter(P.Value <= ALPHA) |>
    arrange(match(level, names(LEVEL_LABELS)), match(phenotype, PHENOTYPES), P.Value)
  page_key <- paste(hits$level, hits$phenotype)
  link_pages <- split(hits, factor(page_key, unique(page_key))) |>
    unname() |>
    map(\(d) {
      level <- d$level[1]
      p <- d$phenotype[1]
      strips <- sprintf(
        "%s%s\nr %.2f (VEH %.2f, PAS %.2f)\np %.2g, FDR %.2g", short(d$label, 32),
        if_else(level == "protein" & d$universe == "mito", " (mito)", ""),
        d$r, d$r_veh, d$r_pas, d$P.Value, d$adj.P.Val
      )
      long <- pmap(list(d$feature, d$universe, strips), \(f, u, s) {
        tibble(
          strip = s, sample_id = mice, value = feature_matrix(level, u)[f, mice],
          phenotype_value = phenotype_of(p)
        )
      }) |>
        list_rbind() |>
        left_join(select(targets, sample_id, group), by = "sample_id") |>
        mutate(group = factor(group, FIGURE_GROUPS), strip = strip_factor(strip))
      effect <- filter(effects, phenotype == p)
      pages <- link_grid(long,
        title = str_glue(
          "{label}: {LEVEL_LABELS[level]} associated with {PHENOTYPE_LABELS[p]}, ",
          "every association at p ≤ 0.05 (n = {nrow(d)})"
        ),
        subtitle = sprintf(
          paste0(
            "Phenotype effect of PAS in this sex: %+.0f%% (p %.2g)%s.\n",
            "Lines within vehicle and within PAS (5 mice each); dashed: all 10 mice, ",
            "the tested association"
          ),
          effect$pct_change, effect$P.Value,
          if (effect$P.Value > ALPHA) {
            ", so these associations follow variation between mice, not a PAS response"
          } else {
            ""
          }
        )
      )
      set_names(
        pages,
        sprintf(
          "association grid, %s, %s, page %d", LEVEL_LABELS[level], PHENOTYPE_LABELS[p],
          seq_along(pages)
        )
      )
    }) |>
    unlist(recursive = FALSE)
  plot_index <- write_report(file.path(out$reports, "03_associations.pdf"), c(
    list(
      "phenotype effects" = effects_page, "phenotype association counts by level" = grid_page,
      "mito-share check" = check_page
    ),
    set_names(r_heatmaps, paste("r heatmap,", names(r_heatmaps))),
    set_names(examples_page, "robust per-phenotype examples"),
    link_pages
  ))

  saveRDS(
    list(
      effects = effects, associations = links, association_summary = link_summary, fgsea = fgsea,
      examples = examples, treatment_r = r_group
    ),
    file.path(out$data, "associations.rds")
  )
  write_workbook(
    list(
      phenotype_effects = effects, association_summary = link_summary,
      protein_associations = filter(links, level == "protein"),
      score_associations = filter(links, level != "protein"),
      pathway_association_fgsea = fgsea, treatment_r = r_group, examples = examples,
      plot_index = plot_index
    ),
    c(
      "Phenotype effects: PAS vs vehicle per phenotype in this sex, % change, 95% CI, Welch p.",
      paste(
        "Per level, universe and phenotype: associations at p ≤ 0.05, chance, FDR tiers;",
        "adjusted rows add pellet mito share."
      ),
      paste(
        "Every protein's association with each phenotype across this sex's ten mice:",
        "slope (log2 per SD), 95% CI, t, p, BH FDR, r pooled, r in VEH, r in PAS."
      ),
      paste(
        "Every pathway, module, complex and cluster score's association with each phenotype,",
        "same columns."
      ),
      paste(
        "fgsea on protein association t per phenotype (gene permutation),",
        "beside the same set's singscore association (mouse level)."
      ),
      paste(
        "Treatment effect beside the associations: every feature's r with PAS (PAS = 1)",
        "in this sex, with the contrast's p and FDR."
      ),
      "The robust example protein per phenotype.",
      "Which PDF page shows what."
    ),
    out$workbook, inputs
  )
  link_summary
}

# ROC curve points for one feature, drawn on the side it separates.
roc_points <- function(x, case, auc) {
  r <- pROC::roc(
    factor(case, c(FALSE, TRUE)), x,
    direction = if (auc >= 0.5) "<" else ">", quiet = TRUE
  )
  tibble(fpr = 1 - r$specificities, tpr = r$sensitivities) |> arrange(fpr, tpr)
}

CLASSIFICATION_LEVELS <- c(
  pathway = "Pathways", module = "Modules", complex = "Complexes", cluster = "Clusters"
)

# Classification: which single pathway scores, modules, complexes and clusters separate PAS from
# vehicle mice in this sex, sliced from 09_Classification (AUC, DeLong CI, exact Wilcoxon p, BH
# within level, perfect separations against all 252 relabellings). Female and male arcs only.
arc_classification <- function(arc) {
  sex <- ARCS[[arc]]$sex
  stopifnot(!is.na(sex))
  out <- arc_paths(arc, "04_Classification")
  inputs <- c(classification = trunk_path("09_Classification", "classification.rds"))
  cls <- readRDS(inputs[["classification"]])
  label <- ARCS[[arc]]$label
  targets <- readRDS(trunk_path("07_Fit", "fit.rds"))$targets
  keep <- targets$sex == sex
  case <- targets$treatment[keep] == "PAS"
  roc <- filter(cls$roc, arc == !!arc) |>
    select(-arc)
  perfect <- filter(cls$perfect, arc == !!arc) |>
    select(-arc)
  summary <- filter(cls$summary, arc == !!arc) |>
    select(-arc)

  perfect_page <- perfect |>
    mutate(level = factor(CLASSIFICATION_LEVELS[level], rev(CLASSIFICATION_LEVELS))) |>
    ggplot(aes(observed, level)) +
    geom_errorbar(
      aes(xmin = 0, xmax = null_q95),
      width = 0.3, colour = "grey60", orientation = "y"
    ) +
    geom_point(aes(x = null_median), shape = 124, size = 3, colour = "grey40") +
    geom_point(size = 2, colour = arc_colour(arc)) +
    geom_text(aes(label = sprintf("p %.3f", p)), nudge_y = 0.3, hjust = 0, size = 2) +
    labs(
      x = "Features separating PAS from vehicle perfectly", y = NULL,
      title = str_glue("{label}: perfect separations against chance"),
      subtitle = paste(
        "Dot: observed; grey: median and 95th percentile",
        "over all 252 relabellings of the same mice"
      )
    ) +
    theme_figure()
  auc_page <- roc |>
    mutate(level = factor(CLASSIFICATION_LEVELS[level], CLASSIFICATION_LEVELS)) |>
    ggplot(aes(auc)) +
    geom_histogram(breaks = seq(0, 1, 0.05), fill = "grey65", colour = "white", linewidth = 0.2) +
    facet_wrap(~level, scales = "free_y", ncol = 1) +
    labs(
      x = "AUC (1: higher in PAS)", y = "Features",
      title = str_glue("{label}: AUC of every feature")
    ) +
    theme_figure()
  curve_pages <- roc |>
    filter(p <= ALPHA) |>
    group_by(level) |>
    group_split() |>
    map(\(d) {
      level <- d$level[1]
      m <- cls$features[[level]][, keep, drop = FALSE]
      long <- pmap(d, \(feature, auc, p, fdr, r_group, ...) {
        roc_points(m[feature, ], case, auc) |>
          mutate(strip = sprintf(
            "%s\nAUC %.2f (%s), r %.2f\np %.3f, FDR %.2g", short(feature, 24), auc,
            if (auc >= 0.5) "higher" else "lower", r_group, p, fdr
          ))
      }) |>
        list_rbind() |>
        mutate(strip = strip_factor(strip))
      pages <- paginate(long, \(x) {
        ggplot(x, aes(fpr, tpr)) +
          geom_abline(linetype = "dashed", linewidth = 0.3, colour = "grey60") +
          geom_step(direction = "vh", colour = arc_colour(arc), linewidth = 0.5) +
          facet_wrap(~strip, ncol = 5, nrow = 4, drop = FALSE) +
          coord_equal() +
          labs(
            x = "1 - specificity", y = "Sensitivity",
            title = str_glue(
              "{label}: {CLASSIFICATION_LEVELS[level]}, PAS vs vehicle, ",
              "every feature at p ≤ 0.05 (n = {nrow(d)})"
            ),
            subtitle = paste(
              "Exact Wilcoxon p, BH within level; higher or lower in PAS;",
              "r with PAS (PAS = 1)"
            )
          ) +
          theme_figure() +
          theme(
            strip.text = element_text(size = 4.6, face = "plain"),
            axis.text = element_text(size = 4.5)
          )
      })
      set_names(
        pages,
        sprintf("ROC grid, %s, page %d", CLASSIFICATION_LEVELS[level], seq_along(pages))
      )
    }) |>
    unlist(recursive = FALSE)
  plot_index <- write_report(file.path(out$reports, "04_classification.pdf"), c(
    list("perfect separations against chance" = perfect_page, "AUC distributions" = auc_page),
    curve_pages
  ))

  saveRDS(
    list(roc = roc, perfect = perfect, summary = summary),
    file.path(out$data, "classification.rds")
  )
  write_workbook(
    list(summary = summary, roc = roc, perfect_separations = perfect, plot_index = plot_index),
    c(
      "Per level: features, at p ≤ 0.05, expected by chance, FDR ≤ 0.05, perfect separations.",
      paste(
        "Every feature, by p: AUC (PAS vs vehicle), DeLong 95% CI (empty at perfect separation),",
        "exact Wilcoxon p, BH FDR within level, r with PAS."
      ),
      "Per level: perfect separations against all 252 relabellings.",
      "Which PDF page shows what."
    ),
    out$workbook, inputs
  )
  summary
}
