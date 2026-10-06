# The steps every arc runs on its own slice of the trunk, so the female, male and sex-difference
# arcs differ only in the contrast and mice they read. Each step writes one workbook, one report
# PDF and the rds its figures read. Sourced by the arc scripts; it sources R/helpers.R.

source(here::here("R", "helpers.R"))

trunk_path <- \(step, file) here("00_Trunk", step, "c_data", file)

arc_paths <- function(arc, step) {
  out <- stage_paths(ARCS[[arc]]$dir, step)
  out$workbook <- file.path(out$data, paste0(tolower(step), ".xlsx"))
  out
}

# The mice an arc looks at: one sex for the PAS arcs, all twenty for the sex difference.
arc_mice <- function(arc, targets) {
  sex <- ARCS[[arc]]$sex
  if (is.na(sex)) targets$sample_id else targets$sample_id[targets$sex == sex]
}

arc_colour <- \(arc) unname(CONTRAST_COLOURS[CONTRAST_LABELS[ARCS[[arc]]$contrast]])

tier_of <- function(p, fdr, pi) {
  factor(
    case_when(
      fdr <= ALPHA ~ "FDR ≤ 0.05", p <= ALPHA & pi <= ALPHA ~ "p and Π ≤ 0.05",
      p <= ALPHA ~ "p ≤ 0.05", .default = "p > 0.05"
    ),
    c("p > 0.05", "p ≤ 0.05", "p and Π ≤ 0.05", "FDR ≤ 0.05")
  )
}
tier_colours <- function(arc) {
  colour <- arc_colour(arc)
  c(
    "p > 0.05" = "grey82", "p ≤ 0.05" = colorspace::lighten(colour, 0.55),
    "p and Π ≤ 0.05" = colour, "FDR ≤ 0.05" = "black"
  )
}

# One mouse per point, grouped and coloured by group, for a handful of features.
per_mouse_features <- function(m, targets, features, labels, title, subtitle) {
  as_tibble(m[features, , drop = FALSE], rownames = "feature") |>
    pivot_longer(-feature, names_to = "sample_id") |>
    left_join(select(targets, sample_id, group, sex), by = "sample_id") |>
    mutate(
      group = factor(group, FIGURE_GROUPS),
      feature = factor(labels[feature], unique(labels[features]))
    ) |>
    ggplot(aes(group, value, colour = group, shape = sex)) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.25, colour = "grey40") +
    geom_point(position = position_jitter(width = 0.12, seed = SEED), size = 1.3) +
    facet_wrap(~feature, scales = "free_y") +
    scale_colour_manual(values = GROUP_COLOURS, guide = "none") +
    scale_shape_manual(values = SEX_SHAPES, guide = "none") +
    labs(x = NULL, y = "log2 abundance", title = title, subtitle = subtitle) +
    theme_figure()
}

# Proteins: the arc's contrast sliced from the shared fit, its tiers, the lists at p and Π, and
# whether those lists are more connected in STRING than chance (descriptive: proteins of one
# complex vary together, so connected lists do not show the effects are real).
arc_proteins <- function(arc) {
  out <- arc_paths(arc, "01_Proteins")
  inputs <- c(fit = trunk_path("07_Fit", "fit.rds"))
  fit <- readRDS(inputs[["fit"]])
  contrast <- ARCS[[arc]]$contrast
  mice <- arc_mice(arc, fit$targets)

  results <- fit$protein_contrasts |>
    filter(contrast == !!contrast) |>
    mutate(tier = tier_of(P.Value, adj.P.Val, pi_score), direction = if_else(logFC > 0, "up", "down"))
  summary <- results |>
    reframe(tier_counts(P.Value, adj.P.Val, pi_score), .by = c(universe, direction))
  lists <- results |>
    filter(universe == "full", P.Value <= ALPHA, pi_score <= ALPHA, !is.na(gene), gene != "") |>
    select(direction, feature, gene, logFC, P.Value, adj.P.Val, pi_score) |>
    arrange(direction, pi_score)

  string_db <- STRINGdb::STRINGdb$new(
    version = "12.0", species = 10090, score_threshold = 400, network_type = "full",
    input_directory = here("00_Trunk", "04_Gene_Sets", "cache")
  )
  detected <- unique(filter(results, universe == "full", !is.na(gene), gene != "")$gene)
  string_map <- suppressWarnings(string_db$map(data.frame(gene = detected), "gene",
    removeUnmappedRows = TRUE
  ))
  string_db$set_background(unique(string_map$STRING_id))
  string <- c(split(lists$gene, lists$direction),
    `FDR ≤ 0.05` = list(filter(results, universe == "full", adj.P.Val <= ALPHA)$gene)
  ) |>
    imap(\(genes, list_name) {
      ids <- unique(string_map$STRING_id[string_map$gene %in% genes])
      e <- if (length(ids) >= 2) string_db$get_ppi_enrichment(ids) else list(edges = 0, lambda = 0, enrichment = 1)
      tibble(
        list = list_name, proteins = length(genes), mapped = length(ids), edges = e$edges,
        expected = e$lambda, fold = if_else(e$lambda > 0, e$edges / e$lambda, NA_real_),
        p = e$enrichment
      )
    }) |>
    list_rbind() |>
    mutate(fdr = p.adjust(p, "BH"))

  label <- ARCS[[arc]]$label
  volcano <- results |>
    mutate(universe = if_else(universe == "full", "Full proteome", "Mitochondrial universe")) |>
    ggplot(aes(logFC, -log10(P.Value), colour = tier)) +
    geom_point(size = 0.6) +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", linewidth = 0.3) +
    ggrepel::geom_text_repel(
      data = \(d) slice_min(filter(d, P.Value <= ALPHA), pi_score, n = 10, by = universe),
      aes(label = gene), size = 2, colour = "black", max.overlaps = Inf, min.segment.length = 0,
      seed = SEED
    ) +
    facet_wrap(~universe) +
    scale_colour_manual(values = tier_colours(arc), name = NULL, drop = FALSE) +
    labs(
      x = "log2 fold change", y = "-log10 p", title = str_glue("{label}: {CONTRAST_LABELS[contrast]}"),
      subtitle = "limpa dpcDE, BH within universe; labels: 10 smallest Π among p ≤ 0.05"
    ) +
    theme_figure() +
    theme(legend.position = "bottom")
  counts <- summary |>
    filter(universe == "full") |>
    pivot_longer(c(p_le_0.05, p_and_pi_le_0.05, fdr_le_0.05), names_to = "tier", values_to = "n") |>
    mutate(
      tier = factor(tier, c("p_le_0.05", "p_and_pi_le_0.05", "fdr_le_0.05"), c("p ≤ 0.05", "p and Π ≤ 0.05", "FDR ≤ 0.05")),
      n = if_else(direction == "down", -n, n)
    ) |>
    ggplot(aes(n, tier, fill = direction)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = c(-1, 1) * ALPHA * nrow(fit$E_norm$full) / 2, linetype = "dashed", linewidth = 0.3) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_text(aes(label = abs(n), hjust = if_else(n < 0, 1.2, -0.2)), size = 2.2) +
    scale_fill_manual(values = c(down = "#2166AC", up = "#B2182B"), name = NULL) +
    scale_x_continuous(labels = abs, expand = expansion(mult = 0.15)) +
    labs(
      x = "Proteins (down | up)", y = NULL, title = "Proteins per tier",
      subtitle = sprintf("Full proteome, %d proteins; dashed: p ≤ 0.05 by chance per side", nrow(fit$E_norm$full))
    ) +
    theme_figure()
  hits <- filter(results, universe == "full", P.Value <= ALPHA, pi_score <= ALPHA)
  heatmap <- group_heatmap(fit$E_norm$full[hits$feature, mice, drop = FALSE],
    filter(fit$targets, sample_id %in% mice),
    title = str_glue("{nrow(hits)} proteins at p and Π ≤ 0.05"),
    subtitle = "Rows z-scored across these mice, clustered on 1 - Pearson r (display only)"
  )
  top <- slice_min(hits, pi_score, n = 12, with_ties = FALSE)
  dots <- per_mouse_features(fit$E_norm$full[, mice], filter(fit$targets, sample_id %in% mice),
    top$feature, set_names(top$gene, top$feature),
    title = "Twelve smallest Π among p ≤ 0.05", subtitle = "log2 abundance per mouse; bar: group median"
  )
  plot_index <- write_report(file.path(out$reports, "01_proteins.pdf"), list(
    "volcanoes, full and mitochondrial universe" = volcano, "proteins per tier" = counts,
    "heatmap of p and Π proteins" = heatmap, "per-mouse values, top 12" = dots
  ))

  saveRDS(list(results = results, summary = summary, lists = lists, string = string, mice = mice),
    file.path(out$data, "proteins.rds")
  )
  write_workbook(
    list(
      summary = summary, results = select(results, -genes), lists = lists, string = string,
      plot_index = plot_index
    ),
    c(
      "Per universe and direction: p ≤ 0.05, chance, p and Π ≤ 0.05, FDR tiers.",
      "Every protein in both universes: log2 FC with 95% CI, t, p, BH FDR, Π and tier.",
      "Full-universe proteins at p and Π ≤ 0.05, by direction.",
      "STRING v12 connectivity of each list against the detected proteome (descriptive).",
      "Which PDF page shows what."
    ),
    out$workbook, inputs
  )
  summary
}

# One contrast of a trunk study, so enrichVolcano tests only what the arc asks.
arc_study <- function(arc, universe) {
  study <- suppressMessages(enrichVolcano::read_study(
    trunk_path("07_Fit", paste0("study_", universe)),
    species = "Mus musculus"
  ))
  contrast <- ARCS[[arc]]$contrast
  study$da <- study$da[study$da$contrast == contrast, ]
  study$contrasts <- study$contrasts[study$contrasts$name == contrast, , drop = FALSE]
  study
}

flat_results <- function(enrichment, universe) {
  imap(enrichment, \(e, test) mutate(e@results, test = test, universe = universe, .before = 1)) |>
    list_rbind() |>
    mutate(leading_edge = map_chr(leading_edge, \(x) paste(x, collapse = ";")))
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

# Pathways: fgsea (NES, the lead result) and fry on the arc's contrast, GO:BP in the full proteome
# and MitoCarta in the mitochondrial universe, BH within each collection over every term tested.
# GO Slim generic BP families (every slim ancestor a term reaches) organise the GO:BP terms; they
# never choose what is tested. Rings take the 12 strongest terms per universe after the
# EnrichmentMap collapse; the tree takes the 20 strongest uncollapsed GO:BP terms. Display only.
arc_pathways <- function(arc) {
  out <- arc_paths(arc, "02_Pathways")
  inputs <- c(
    gene_sets = trunk_path("04_Gene_Sets", "gene_sets.rds"),
    study_full = trunk_path("07_Fit", "study_full/da_results.csv"),
    study_mito = trunk_path("07_Fit", "study_mito/da_results.csv")
  )
  gs <- readRDS(inputs[["gene_sets"]])
  contrast <- ARCS[[arc]]$contrast
  collections <- list(full = gs$full["GO:BP"], mito = gs$mito)
  min_size <- c(full = 15, mito = 10)
  studies <- map(set_names(names(collections)), \(u) arc_study(arc, u))

  enrichment <- imap(collections, \(sets, u) {
    set.seed(SEED)
    en <- suppressMessages(enrichVolcano::run_enrichment(studies[[u]], sets,
      tests = c("fgsea", "fry"), min_size = min_size[[u]], max_size = 500
    ))
    map(en, \(e) enrichVolcano::dedup_terms(e, unlist(unname(sets), recursive = FALSE)))
  })
  set.seed(SEED)
  compartments <- suppressMessages(enrichVolcano::run_enrichment(studies$full, gs$compartments,
    tests = c("fgsea", "fry"), min_size = 15, max_size = Inf
  ))

  families <- gs$goslim |>
    transmute(term, family = str_split(slim, "; ")) |>
    unnest(family) |>
    filter(family != "")
  results <- imap(enrichment, flat_results) |>
    list_rbind() |>
    left_join(summarise(families, families = paste(family, collapse = "; "), .by = term), by = "term") |>
    select(
      test, universe, collection = database, term, NES = score, p, padj, size, direction,
      dedup_status, families, leading_edge
    )
  compartment_results <- imap(compartments, \(e, test) mutate(e@results, test = test)) |>
    list_rbind() |>
    transmute(test, compartment = term, NES = score, p, padj, size, direction)
  summary <- results |>
    reframe(tier_counts(p, padj), .by = c(test, universe, collection)) |>
    select(-p_and_pi_le_0.05)

  fgsea <- filter(results, test == "fgsea")
  family_summary <- fgsea |>
    filter(collection == "GO:BP") |>
    inner_join(families, by = "term", relationship = "many-to-many") |>
    summarise(
      terms = n(), up_p = sum(p <= ALPHA & NES > 0), down_p = sum(p <= ALPHA & NES < 0),
      up_fdr = sum(padj <= ALPHA & NES > 0), down_fdr = sum(padj <= ALPHA & NES < 0),
      .by = family
    ) |>
    arrange(desc(up_p + down_p))

  display <- fgsea |>
    filter(!(collection == "GO:BP" & term %in% gs$slim_terms))
  ring_terms <- display |>
    filter(dedup_status == "kept" | is.na(dedup_status)) |>
    slice_min(p, n = 12, by = universe, with_ties = FALSE)
  tree_terms <- display |>
    filter(collection == "GO:BP") |>
    slice_min(p, n = 20, with_ties = FALSE)
  sets <- gs$full[["GO:BP"]][tree_terms$term]
  tree <- hclust(as.dist(1 - enrichmentmap_similarity(sets)), method = "average")
  tree_terms <- tree_terms |>
    mutate(cluster = cutree(tree, h = 1 - 0.375)[term]) |>
    mutate(cluster_name = term[which.min(p)], .by = cluster)
  literature <- fgsea |>
    inner_join(gs$literature, by = c(term = "set")) |>
    mutate(list_fdr = p.adjust(p, "BH")) |>
    select(hypothesis, expected, collection, term, NES, p, padj, list_fdr)

  label <- ARCS[[arc]]$label
  rings <- imap(c(full = "GO:BP", mito = "MitoCarta"), \(db, u) {
    terms <- filter(ring_terms, universe == u)
    enrichVolcano::plot_volcano_ring(studies[[u]]$da, enrichment[[u]]$fgsea,
      contrast = contrast, databases = db, terms = terms$term,
      labels = set_names(
        paste0(enrichVolcano::clean_label(terms$term, width = 22), if_else(terms$padj <= ALPHA, " *", "")),
        terms$term
      ),
      title = str_glue("{label}: {db}"),
      subtitle = "12 strongest terms after the EnrichmentMap collapse; * FDR ≤ 0.05",
      p_threshold = ALPHA, show_counts = FALSE, score_limits = c(-3, 3), arc_span = "auto",
      theme = enrichVolcano::plot_theme(base_size = 7, base_family = "Helvetica")
    )
  })
  family_page <- family_summary |>
    slice_max(up_p + down_p, n = 20, with_ties = FALSE) |>
    pivot_longer(c(up_p, down_p), names_to = "direction", values_to = "n") |>
    mutate(
      n = if_else(direction == "down_p", -n, n),
      family = forcats::fct_reorder(str_to_sentence(family), abs(n), sum)
    ) |>
    ggplot(aes(n, family, fill = direction)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = c(down_p = "#2166AC", up_p = "#B2182B"), labels = c("down", "up"), name = NULL) +
    scale_x_continuous(labels = abs) +
    labs(
      x = "GO:BP terms at p ≤ 0.05 (down | up)", y = NULL, title = "Signal by GO Slim family",
      subtitle = "fgsea; a term counts in every family it belongs to; 20 families with most terms"
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
    scale_fill_gradient2(low = "#2166AC", high = "#B2182B", limits = c(-3, 3), oob = scales::squish) +
    labs(
      x = "NES", y = NULL, title = "20 strongest GO:BP terms, ordered by shared genes",
      subtitle = "EnrichmentMap similarity, average linkage; black ring FDR ≤ 0.05"
    ) +
    theme_figure()
  compartment_page <- compartment_results |>
    filter(test == "fry") |>
    mutate(
      signed = -log10(p) * if_else(direction == "up", 1, -1),
      compartment = forcats::fct_reorder(enrichVolcano::clean_label(compartment, width = 40), signed)
    ) |>
    ggplot(aes(signed, compartment, fill = padj <= ALPHA)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = c(-1, 1) * -log10(ALPHA), linetype = "dashed", linewidth = 0.3) +
    scale_fill_manual(values = c(`TRUE` = arc_colour(arc), `FALSE` = "grey70"), labels = c("FDR > 0.05", "FDR ≤ 0.05"), name = NULL) +
    labs(
      x = "Signed -log10 p (fry)", y = NULL, title = "GO Slim cellular compartments",
      subtitle = "fry rotation test, full proteome; dashed p = 0.05"
    ) +
    theme_figure()
  plot_index <- write_report(file.path(out$reports, "02_pathways.pdf"), list(
    "ring, GO:BP" = rings$full, "ring, MitoCarta" = rings$mito, "signal by GO Slim family" = family_page,
    "strongest GO:BP terms by shared genes" = tree_page, "compartments (fry)" = compartment_page
  ))

  saveRDS(
    list(
      enrichment = enrichment, compartments = compartment_results, results = results,
      ring_terms = ring_terms, tree = tree, tree_terms = tree_terms, families = family_summary,
      literature = literature
    ),
    file.path(out$data, "pathways.rds")
  )
  write_workbook(
    list(
      summary = summary, fgsea = filter(results, test == "fgsea"), fry = filter(results, test == "fry"),
      families = family_summary, compartments = compartment_results,
      ring_terms = select(ring_terms, universe, collection, term, NES, p, padj),
      tree_terms = select(tree_terms, term, NES, p, padj, cluster, cluster_name),
      literature = literature, plot_index = plot_index
    ),
    c(
      "Per test, universe and collection: p ≤ 0.05, chance, FDR tiers.",
      "fgsea per term: NES, p, BH FDR within collection, size, collapse status, GO Slim families, leading edge.",
      "fry per term: direction, p, BH FDR within collection.",
      "Per GO Slim family: GO:BP terms tested and at p and FDR ≤ 0.05 by direction (fgsea).",
      "GO Slim cellular compartments: fgsea and fry, full proteome.",
      "Ring terms: 12 strongest per universe after the EnrichmentMap collapse.",
      "Tree terms: 20 strongest GO:BP terms, their cluster at EnrichmentMap 0.375 and its lead term.",
      "Literature-informed sets present in GO:BP or MitoCarta, with BH across the list.",
      "Which PDF page shows what."
    ),
    out$workbook, inputs
  )
  summary
}

# The protein with the smallest link p whose slope keeps its sign and p <= 0.05 when any one mouse
# is left out (ordinary least squares within the sex), so an example is not one mouse's doing.
robust_example <- function(links, E, z, mice, max_tries = 50) {
  candidates <- slice_min(links, P.Value, n = max_tries, with_ties = FALSE)
  for (i in seq_len(nrow(candidates))) {
    y <- E[candidates$feature[i], mice]
    holds <- map_lgl(seq_along(mice), \(k) {
      f <- summary(lm(y[-k] ~ z[-k]))$coefficients
      sign(f[2, 1]) == sign(candidates$estimate[i]) && f[2, 4] <= ALPHA
    })
    if (all(holds)) {
      return(candidates[i, ])
    }
  }
  candidates[0, ]
}

# Associations: what PAS did to each phenotype in this sex, and which proteins, GO:BP and MitoCarta
# scores, modules, complexes and clusters move with each phenotype across this sex's PAS and
# vehicle mice (07_Fit). Pathway patterns come from fgsea on the protein link t (NES, a gene-
# permutation test: exploratory); the mouse-level test of the same term is its singscore link.
arc_associations <- function(arc) {
  out <- arc_paths(arc, "03_Associations")
  inputs <- c(
    fit = trunk_path("07_Fit", "fit.rds"),
    gene_sets = trunk_path("04_Gene_Sets", "gene_sets.rds")
  )
  fit <- readRDS(inputs[["fit"]])
  gs <- readRDS(inputs[["gene_sets"]])
  sex <- ARCS[[arc]]$sex
  label <- ARCS[[arc]]$label
  mice <- arc_mice(arc, fit$targets)
  targets <- filter(fit$targets, sample_id %in% mice)

  effects <- filter(fit$phenotype_effects, sex == !!sex)
  protein_links <- filter(fit$protein_links, sex == arc)
  score_links <- filter(fit$score_links, sex == arc)
  links <- bind_rows(mutate(protein_links, level = "protein"), score_links) |>
    filter(!adjusted)
  link_summary <- bind_rows(mutate(protein_links, level = "protein"), score_links) |>
    reframe(tier_counts(P.Value, adj.P.Val), .by = c(level, universe, phenotype, adjusted)) |>
    select(-p_and_pi_le_0.05)

  collections <- list(full = gs$full["GO:BP"], mito = gs$mito)
  min_size <- c(full = 15, mito = 10)
  fgsea <- imap(collections, \(sets, u) {
    da <- protein_links |>
      filter(!adjusted, universe == u) |>
      transmute(
        protein = feature, gene, contrast = phenotype, logFC = estimate, t, P.Value, adj.P.Val,
        abundance = AveExpr
      ) |>
      enrichVolcano::as_da(species = "Mus musculus")
    set.seed(SEED)
    suppressMessages(enrichVolcano::run_enrichment(da, sets,
      tests = "fgsea", min_size = min_size[[u]], max_size = 500
    ))$fgsea@results |>
      mutate(universe = u, .before = 1)
  }) |>
    list_rbind() |>
    transmute(universe, collection = database, phenotype = contrast, term, NES = score, p, padj,
      size, leading_edge = map_chr(leading_edge, \(x) paste(x, collapse = ";"))
    )
  matched <- fgsea |>
    left_join(
      score_links |>
        filter(!adjusted, level == "pathway") |>
        select(phenotype, term = feature, score_slope = estimate, score_t = t, score_p = P.Value, score_fdr = adj.P.Val),
      by = c("phenotype", "term")
    )

  examples <- map(PHENOTYPES, \(p) {
    z <- as.numeric(scale(fit$phenotypes[[p]]))[match(mice, fit$targets$sample_id)]
    protein_links |>
      filter(!adjusted, universe == "full", phenotype == p) |>
      robust_example(fit$E_norm$full, z, mice)
  }) |>
    list_rbind()

  effects_page <- phenotype_effect_plot(effects,
    title = str_glue("{label}: PAS effect on each phenotype"),
    subtitle = "% change against vehicle, Welch t on log values, n = 5 v 5"
  )
  grid_page <- link_summary |>
    filter(!adjusted) |>
    mutate(
      level = factor(paste(level, if_else(universe == "mito", "(mito)", "")), c(
        "protein ", "protein (mito)", "pathway ", "pathway (mito)", "module ", "complex ", "cluster "
      ), c("Proteins", "Proteins, mito", "GO:BP", "MitoCarta", "Modules", "Complexes", "Clusters")),
      trait = factor(PHENOTYPE_LABELS[phenotype], rev(PHENOTYPE_LABELS)),
      ratio = p_le_0.05 / expected_by_chance
    ) |>
    ggplot(aes(level, trait, fill = ratio)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = str_glue("{p_le_0.05}\nFDR {fdr_le_0.05}")), size = 2) +
    scale_fill_gradient(low = "white", high = "grey35", name = "p ≤ 0.05 /\nchance") +
    labs(
      x = NULL, y = NULL, title = str_glue("{label}: features moving with each phenotype"),
      subtitle = "Across PAS and vehicle mice of this sex (n = 10); BH within level and phenotype"
    ) +
    theme_figure()
  top_terms <- fgsea |>
    filter(collection == "GO:BP", !term %in% gs$slim_terms) |>
    summarise(best = min(p), .by = term) |>
    slice_min(best, n = 30, with_ties = FALSE) |>
    pull(term)
  nes <- matched |>
    filter(term %in% top_terms) |>
    select(term, phenotype, NES) |>
    pivot_wider(names_from = phenotype, values_from = NES)
  order <- nes$term[hclust(dist(as.matrix(nes[-1])), "average")$order]
  heatmap_page <- matched |>
    filter(term %in% top_terms) |>
    mutate(
      term = factor(term, order),
      trait = factor(PHENOTYPE_LABELS[phenotype], PHENOTYPE_LABELS)
    ) |>
    ggplot(aes(trait, term, fill = NES)) +
    geom_tile(colour = "white") +
    geom_point(data = \(d) filter(d, score_fdr <= ALPHA), size = 1) +
    scale_fill_gradient2(low = "#2166AC", high = "#B2182B", limits = c(-3, 3), oob = scales::squish) +
    scale_y_discrete(labels = \(x) enrichVolcano::clean_label(x, width = 60)) +
    labs(
      x = NULL, y = NULL, title = str_glue("{label}: GO:BP terms moving with each phenotype"),
      subtitle = "Colour: fgsea NES on protein link t (exploratory); dot: singscore link FDR ≤ 0.05; rows clustered"
    ) +
    theme_figure() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  example_page <- examples |>
    rowwise() |>
    reframe(
      phenotype, gene, P.Value, adj.P.Val,
      sample_id = mice, abundance = fit$E_norm$full[feature, mice],
      value = fit$phenotypes[[phenotype]][match(mice, fit$targets$sample_id)]
    ) |>
    left_join(select(targets, sample_id, group), by = "sample_id") |>
    mutate(strip = factor(
      sprintf("%s · %s\np %.2g, FDR %.2g", PHENOTYPE_LABELS[phenotype], gene, P.Value, adj.P.Val),
      unique(sprintf("%s · %s\np %.2g, FDR %.2g", PHENOTYPE_LABELS[phenotype], gene, P.Value, adj.P.Val))
    )) |>
    ggplot(aes(value, abundance)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = arc_colour(arc), linewidth = 0.5) +
    geom_point(aes(colour = factor(group, FIGURE_GROUPS)), size = 1.6) +
    facet_wrap(~strip, scales = "free", nrow = 2) +
    scale_colour_manual(values = GROUP_COLOURS, name = NULL) +
    labs(
      x = "Phenotype", y = "log2 abundance",
      title = str_glue("{label}: strongest robust protein per phenotype"),
      subtitle = "Smallest link p that keeps its sign and p ≤ 0.05 with any one mouse left out"
    ) +
    theme_figure()
  plot_index <- write_report(file.path(out$reports, "03_associations.pdf"), list(
    "phenotype effects" = effects_page, "link counts by level" = grid_page,
    "GO:BP link heatmap" = heatmap_page, "robust per-mouse examples" = example_page
  ))

  saveRDS(
    list(effects = effects, links = links, link_summary = link_summary, fgsea = matched, examples = examples),
    file.path(out$data, "associations.rds")
  )
  write_workbook(
    list(
      phenotype_effects = effects, link_summary = link_summary,
      protein_links = filter(links, level == "protein"), score_links = filter(links, level != "protein"),
      pathway_links = matched, examples = examples, plot_index = plot_index
    ),
    c(
      "PAS vs vehicle per phenotype in this sex: % change, 95% CI, Welch p.",
      "Per level, universe and phenotype: p ≤ 0.05, chance, FDR tiers; adjusted rows add pellet mito share.",
      "Every protein's link with each phenotype (log2 per SD of the phenotype), unadjusted.",
      "Every GO:BP, MitoCarta, module, complex and cluster score's link with each phenotype.",
      "fgsea on protein link t per phenotype, matched to the same term's singscore link.",
      "The robust example protein per phenotype.",
      "Which PDF page shows what."
    ),
    out$workbook, inputs
  )
  link_summary
}

# Every feature an arc classifies with, one matrix per level, rows labelled for display.
arc_features <- function(fit) {
  genes <- distinct(fit$protein_contrasts, feature, gene)
  collection <- distinct(fit$score_contrasts, feature, level)
  level_rows <- \(lvl) collection$feature[collection$level == lvl]
  pathways <- rbind(
    fit$scores$full[intersect(level_rows("pathway"), rownames(fit$scores$full)), ],
    fit$scores$mito
  )
  list(
    protein = fit$E_norm$full, pathway = pathways,
    complex = fit$scores$full[intersect(level_rows("complex"), rownames(fit$scores$full)), ],
    cluster = fit$scores$full[intersect(level_rows("cluster"), rownames(fit$scores$full)), ],
    module = fit$scores$module
  ) |>
    map(\(m) m[, fit$targets$sample_id, drop = FALSE]) |>
    structure(gene = set_names(genes$gene, genes$feature))
}

# How well each feature alone separates two groups of mice: AUC (pROC; 1 = higher in the case
# group) with DeLong's interval, left empty at perfect separation where it is undefined, and the
# exact Wilcoxon p; BH within the level.
roc_table <- function(m, case) {
  y <- factor(case, levels = c(FALSE, TRUE))
  map(rownames(m), \(f) {
    x <- m[f, ]
    roc <- pROC::roc(y, x, levels = levels(y), direction = "<", quiet = TRUE)
    auc <- as.numeric(roc$auc)
    ci <- if (auc %in% c(0, 1)) rep(NA_real_, 3) else as.numeric(pROC::ci.auc(roc, method = "delong"))
    tibble(
      feature = f, auc = auc, auc_low = ci[1], auc_high = ci[3],
      p = wilcox.test(x[case], x[!case], exact = TRUE)$p.value
    )
  }) |>
    list_rbind() |>
    mutate(fdr = p.adjust(p, "BH"))
}

# How many features separate the groups perfectly, against every relabelling of the same mice
# (choose(10, 5) = 252 for five against five): the exact null for that count.
perfect_null <- function(m, case) {
  keys <- \(idx) apply(idx, 1, \(i) paste(sort(i), collapse = ","))
  n <- sum(case)
  ord <- t(apply(m, 1, order))
  top <- keys(ord[, (ncol(m) - n + 1):ncol(m), drop = FALSE])
  bottom <- keys(ord[, seq_len(n), drop = FALSE])
  count <- \(labelled) sum(top == labelled | bottom == labelled)
  null <- map_int(combn(ncol(m), n, simplify = FALSE), \(i) count(paste(sort(i), collapse = ",")))
  observed <- count(paste(which(case), collapse = ","))
  tibble(
    observed = observed, null_median = median(null), null_q95 = unname(quantile(null, 0.95)),
    relabellings = length(null), p = mean(null >= observed)
  )
}

# Sex difference in the PAS response, mouse by mouse: each PAS mouse's value minus the mean of the
# vehicle mice of its own sex, separated by sex (AUC by ranks). The vehicle mean is shared by all
# mice of a sex, so AUC is judged against shuffles of sex within PAS and within vehicle that
# recompute the vehicle means, not against 0.5.
change_roc <- function(m, targets, n_shuffle = 1000) {
  pas <- which(targets$treatment == "PAS")
  veh <- which(targets$treatment == "VEH")
  auc_of <- function(male_pas, male_veh) {
    veh_mean <- cbind(
      rowMeans(m[, veh[!male_veh], drop = FALSE]), rowMeans(m[, veh[male_veh], drop = FALSE])
    )
    change <- m[, pas] - veh_mean[, ifelse(male_pas, 2, 1)]
    ranks <- matrixStats::rowRanks(change)
    n_m <- sum(male_pas)
    (rowSums(ranks[, male_pas, drop = FALSE]) - n_m * (n_m + 1) / 2) / (n_m * sum(!male_pas))
  }
  observed <- auc_of(targets$sex[pas] == "M", targets$sex[veh] == "M")
  set.seed(SEED)
  null <- replicate(n_shuffle, auc_of(
    sample(targets$sex[pas] == "M"), sample(targets$sex[veh] == "M")
  ))
  tibble(
    feature = rownames(m), auc = observed,
    p = (1 + rowSums(abs(null - 0.5) >= abs(observed - 0.5))) / (n_shuffle + 1)
  ) |>
    mutate(fdr = p.adjust(p, "BH"))
}

# Classification: which single proteins, pathway scores, complexes, clusters and modules separate
# the arc's groups. Female and male arcs: PAS against vehicle within sex. Sex difference: male
# against female within vehicle, within PAS, and on each PAS mouse's change from its own sex's
# vehicle mean.
arc_classification <- function(arc) {
  out <- arc_paths(arc, "04_Classification")
  inputs <- c(fit = trunk_path("07_Fit", "fit.rds"))
  fit <- readRDS(inputs[["fit"]])
  features <- arc_features(fit)
  gene <- attr(features, "gene")
  label <- ARCS[[arc]]$label
  targets <- fit$targets
  sex <- ARCS[[arc]]$sex
  groups <- if (is.na(sex)) {
    list(
      `Male vs female, vehicle` = c(case = "VEH_M", control = "VEH_F"),
      `Male vs female, PAS` = c(case = "PAS_M", control = "PAS_F")
    )
  } else {
    set_names(list(c(case = paste0("PAS_", sex), control = paste0("VEH_", sex))), "PAS vs vehicle")
  }

  group_rocs <- imap(groups, \(g, comparison) {
    keep <- targets$group %in% g
    case <- targets$group[keep] == g[["case"]]
    imap(features, \(m, level) {
      list(
        roc = roc_table(m[, keep], case) |> mutate(level = level, comparison = comparison, .before = 1),
        perfect = perfect_null(m[, keep], case) |> mutate(level = level, comparison = comparison, .before = 1)
      )
    })
  }) |>
    unlist(recursive = FALSE)
  roc <- map(group_rocs, "roc") |> list_rbind()
  perfect <- map(group_rocs, "perfect") |> list_rbind()
  if (is.na(sex)) {
    roc <- bind_rows(roc, imap(features, \(m, level) {
      change_roc(m, targets) |> mutate(level = level, comparison = "Male vs female, PAS change", .before = 1)
    }) |> list_rbind())
  }
  roc <- mutate(roc, label = coalesce(gene[feature], feature), .after = feature)

  perfect_page <- perfect |>
    mutate(level = factor(level, names(features))) |>
    ggplot(aes(observed, level)) +
    geom_errorbar(aes(xmin = 0, xmax = null_q95), width = 0.3, colour = "grey60", orientation = "y") +
    geom_point(aes(x = null_median), shape = 124, size = 3, colour = "grey40") +
    geom_point(size = 2, colour = arc_colour(arc)) +
    facet_wrap(~comparison) +
    labs(
      x = "Features separating the groups perfectly", y = NULL,
      title = str_glue("{label}: perfect separations against chance"),
      subtitle = "Dot: observed; grey: median and 95th percentile over all 252 relabellings of the same mice"
    ) +
    theme_figure()
  top <- roc |>
    filter(comparison == first(comparison)) |>
    slice_min(p, n = 2, by = level, with_ties = FALSE)
  keep <- targets$group %in% groups[[1]]
  curve_page <- pmap(top, \(level, feature, label, auc, p, ...) {
    m <- features[[level]]
    case <- targets$group[keep] == groups[[1]][["case"]]
    r <- pROC::roc(factor(case, c(FALSE, TRUE)), m[feature, keep],
      direction = if (auc >= 0.5) "<" else ">", quiet = TRUE
    )
    tibble(
      fpr = 1 - r$specificities, tpr = r$sensitivities,
      strip = sprintf(
        "%s\n%s\nAUC %.2f (%s), p %.3f", str_to_sentence(level),
        str_trunc(enrichVolcano::clean_label(label, width = 200), 28), auc,
        if (auc >= 0.5) "higher" else "lower", p
      )
    ) |>
      arrange(fpr, tpr)
  }) |>
    list_rbind() |>
    ggplot(aes(fpr, tpr)) +
    geom_abline(linetype = "dashed", linewidth = 0.3, colour = "grey60") +
    geom_step(direction = "vh", colour = arc_colour(arc), linewidth = 0.6) +
    facet_wrap(~strip, ncol = 5) +
    coord_equal() +
    labs(
      x = "1 - specificity", y = "Sensitivity",
      title = str_glue("{label}: strongest separators per level, {names(groups)[1]}"),
      subtitle = "Two smallest exact Wilcoxon p per level; higher or lower in the first-named group"
    ) +
    theme_figure()
  auc_page <- roc |>
    mutate(level = factor(level, names(features))) |>
    ggplot(aes(auc)) +
    geom_histogram(breaks = seq(0, 1, 0.05), fill = "grey65", colour = "white", linewidth = 0.2) +
    facet_grid(level ~ comparison, scales = "free_y") +
    labs(x = "AUC (1: higher in the first-named group)", y = "Features", title = str_glue("{label}: AUC of every feature")) +
    theme_figure()
  plot_index <- write_report(file.path(out$reports, "04_classification.pdf"), list(
    "perfect separations against chance" = perfect_page, "ROC curves, strongest separators" = curve_page,
    "AUC distributions" = auc_page
  ))

  saveRDS(list(roc = roc, perfect = perfect, groups = groups), file.path(out$data, "classification.rds"))
  write_workbook(
    list(roc = roc, perfect_separations = perfect, plot_index = plot_index),
    c(
      "Every feature, level and comparison: AUC, DeLong 95% CI (empty at perfect separation), p, BH FDR within level.",
      "Per level and comparison: features separating perfectly, against all 252 relabellings.",
      "Which PDF page shows what."
    ),
    out$workbook, inputs
  )
  count(roc, comparison, level, p_le_0.05 = p <= ALPHA)
}
