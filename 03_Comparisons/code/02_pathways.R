# Step 02 of every comparison: which pathways its contrast shifted.
# Run by 03_Comparisons/<comparison>/a_script/run_comparison.R.

source(here::here("00_Setup", "Shared_functions", "comparisons.R"))

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

# Pathways: the comparison's contrast sliced from 02_Enrichment. fgsea on GO:BP and Reactome in the
# full proteome and MitoCarta in the mitochondrial universe, BH within collection over every tested
# set,
# each beside the same set's singscore contrast; GO Slim compartments and the literature list
# alongside. GO Slim generic BP families (every slim ancestor a term reaches) organise the GO:BP
# terms; they never choose what is tested. Rings take the 12 strongest terms per collection that
# the EnrichmentMap flag keeps; the tree takes the 20 strongest full-proteome terms. Display only.
report_pathways <- function(comparison) {
  out <- comparison_paths(comparison, "02_pathways")
  inputs <- c(
    enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"),
    gene_sets = output_path("01_Processing", "04_Gene_Sets", "gene_sets.rds")
  )
  enr <- readRDS(inputs[["enrichment"]])
  gs <- readRDS(inputs[["gene_sets"]])
  contrast <- COMPARISONS[[comparison]]$contrast
  label <- COMPARISONS[[comparison]]$label
  collections <- c(`GO:BP` = "full", Reactome = "full", MitoCarta = "mito")

  results <- filter(enr$results, contrast == !!contrast) |>
    select(-contrast) |>
    arrange(factor(collection, names(collections)), p)
  tier_table <- filter(enr$summary, contrast == !!contrast) |>
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

  rings <- imap(collections, \(u, db) {
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
        if (is.na(COMPARISONS[[comparison]]$sex)) {
          "\nup: PAS raised the set more in females than in males"
        } else {
          ""
        }
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
      values = c(down_p = DIRECTION_COLOURS[["Down"]], up_p = DIRECTION_COLOURS[["Up"]]),
      labels = c(down_p = "down", up_p = "up"), name = NULL
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
      low = DIRECTION_COLOURS[["Down"]], high = DIRECTION_COLOURS[["Up"]], limits = c(-3, 3),
      oob = scales::squish
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
      values = c(`TRUE` = comparison_colour(comparison), `FALSE` = "grey70"),
      labels = c(`TRUE` = "FDR ≤ 0.05", `FALSE` = "FDR > 0.05"), name = NULL
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
    file.path(out$data, "02_pathways.rds")
  )
  write_workbook(
    list(
      summary = tier_table, fgsea = results, families = family_summary, compartments = compartments,
      ring_terms = select(ring_terms, universe, collection, term, NES, p, padj),
      tree_terms = select(tree_terms, term, NES, p, padj, cluster, cluster_name),
      literature = literature, plot_index = plot_index
    ),
    c(
      "Per universe and collection: p ≤ 0.05, chance, FDR tiers.",
      paste(
        "Every tested set, by p: NES, p, BH FDR within collection, size, redundancy flag,",
        "singscore contrast (log2 FC, p, FDR, r with group), GO Slim families, leading edge."
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
  tier_table
}
