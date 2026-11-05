# Every pathway test the comparisons read, run once: fgsea on each comparison contrast's moderated
# t and on each protein-phenotype association t, camera beside each across mice, side lists, fry
# for one sex's response in the other, ORA of the co-expression clusters in each contrast's Π
# lists, and STRING connectivity. fgsea permutes genes, so camera and the singscore tests from
# 01_Fit are the mouse-level checks on it.

source(here::here("00_Setup", "Shared_functions", "helpers.R"))

out <- stage_paths("02_Statistics", "02_Enrichment")
inputs <- c(
  gene_sets = output_path("01_Processing", "04_Gene_Sets", "gene_sets.rds"),
  design = output_path("01_Processing", "03_Design", "design.rds"),
  fit = output_path("02_Statistics", "01_Fit", "fit.rds"),
  proteins = output_path("01_Processing", "02_Quantification", "proteins.rds")
)
gs <- readRDS(inputs[["gene_sets"]])
fit <- readRDS(inputs[["fit"]])
d <- readRDS(inputs[["design"]])

sets <- map(gs$tested, \(x) unlist(unname(x), recursive = FALSE))
# Each universe's protein contrasts in enrichVolcano's form; it ranks by the moderated t.
da <- map(set_names(names(gs$tested)), \(u) {
  fit$protein_contrasts |>
    filter(universe == u) |>
    transmute(
      protein = feature, gene, contrast, logFC, t, P.Value, adj.P.Val, abundance = AveExpr
    ) |>
    enrichVolcano::as_da(species = "Mus musculus")
})

# One contrast, or one phenotype association, at a time and seeded each time, so a test's result
# does not depend on which others ran beside it.
fgsea_of <- function(da, collection, size) {
  set.seed(SEED)
  suppressMessages(enrichVolcano::run_enrichment(da, collection,
    tests = "fgsea", min_size = size[1], max_size = size[2]
  ))$fgsea@results
}
# The ranking fgsea tested: moderated t, one protein per gene by abundance.
gene_ranks <- function(da) {
  keep <- one_per_gene(da$protein, da$gene, da$abundance)
  set_names(da$rank[keep], da$gene[keep])
}
flatten <- \(x) map_chr(x, \(g) paste(g, collapse = ";"))

ranks <- map(set_names(COMPARISON_CONTRASTS), \(cn) {
  map(da, \(x) gene_ranks(x[x$contrast == cn, ]))
})
# A set's tested members: its genes in the ranking, the same genes in every contrast and
# association.
tested_members <- imap(sets, \(by_set, u) map(by_set, intersect, names(ranks[[1]][[u]])))
# The sets each test shows (the rows sharing `by` are one test): those at FDR <= ALPHA, grouped
# and picked by representative_sets(); the rest, no group and not shown.
mark_shown <- function(tests, by) {
  key <- union(by, c("universe", "term"))
  called <- tests |>
    filter(padj <= ALPHA) |>
    mutate(
      representative_sets(map2(universe, term, \(u, s) tested_members[[u]][[s]]), p),
      .by = all_of(by)
    )
  tests |>
    left_join(select(called, all_of(key), overlap_group, shown), by = key) |>
    mutate(shown = coalesce(shown, FALSE))
}
enrichment_tables <- map(set_names(COMPARISON_CONTRASTS), \(cn) {
  imap(gs$tested, \(collection, u) {
    fgsea_of(da[[u]][da[[u]]$contrast == cn, ], collection, TEST_WINDOW[[u]])
  })
})

families <- gs$goslim |>
  transmute(term, family = str_split(slim, "; ")) |>
  unnest(family) |>
  filter(family != "")
score_stats <- fit$score_contrasts |>
  filter(level == "pathway", contrast %in% COMPARISON_CONTRASTS) |>
  select(
    contrast,
    term = feature, score_logFC = logFC, score_p = P.Value, score_fdr = adj.P.Val, r_group
  )
results <- map(enrichment_tables, \(by_universe) {
  imap(by_universe, \(e, u) mutate(e, universe = u, .before = 1)) |>
    list_rbind()
}) |>
  list_rbind() |>
  mark_shown(c("contrast", "universe")) |>
  left_join(
    summarise(families, families = paste(family, collapse = "; "), .by = term),
    by = "term"
  ) |>
  left_join(score_stats, by = c("contrast", "term")) |>
  transmute(
    contrast, universe,
    collection = database, term, NES = score, p, padj, size, direction,
    overlap_group, shown, score_logFC, score_p, score_fdr, r_group, families,
    leading_edge = flatten(leading_edge)
  ) |>
  arrange(factor(contrast, COMPARISON_CONTRASTS), universe, collection, p)

# The side lists run on fgsea directly, on the ranking the main tests used: enrichVolcano keeps
# fgsea's default nPermSimple, and at that default a set with unbalanced gene-level statistics
# gets no p-value and drops out. fgsea's own advice is nPermSimple = 10000. The windows are
# side_list_size()'s, on the full proteome's ranked genes.
SIDE_LIST_SIZE <- side_list_size(length(ranks[[1]]$full))
side_list <- function(collection, size) {
  side_sets <- unlist(unname(collection), recursive = FALSE)
  map(COMPARISON_CONTRASTS, \(cn) {
    set.seed(SEED)
    fgsea::fgsea(side_sets, ranks[[cn]]$full,
      minSize = size[1], maxSize = size[2], nPermSimple = 10000
    ) |>
      as_tibble() |>
      mutate(contrast = cn)
  }) |>
    list_rbind() |>
    transmute(
      contrast,
      set = pathway, NES, p = pval, padj, size, direction = if_else(NES > 0, "up", "down"),
      leading_edge = flatten(leadingEdge)
    )
}
compartments <- side_list(gs$compartments, SIDE_LIST_SIZE$compartments)
literature <- side_list(list(Literature = gs$literature_sets), SIDE_LIST_SIZE$literature) |>
  rename(list_fdr = padj) |>
  left_join(gs$literature, by = "set") |>
  relocate(hypothesis, expected, .after = contrast)

# Pathway associations: fgsea on each protein's association t with one phenotype within one sex,
# under both association models (genes permuted, so its FDR runs optimistic), beside the same set's
# singscore association under the same model, which is tested across mice. Each phenotype's shown
# sets are picked over both universes together.
association_models <- unique(fit$association_designs$model)
association_fgsea <- expand_grid(
  model = association_models, sex = names(SEX_COMPARISONS), phenotype = PHENOTYPES,
  u = names(gs$tested)
) |>
  pmap(\(model, sex, phenotype, u) {
    fit$protein_associations |>
      filter(model == !!model, sex == !!sex, phenotype == !!phenotype, universe == u) |>
      transmute(
        protein = feature, gene, contrast = phenotype, logFC = estimate, t, P.Value, adj.P.Val,
        abundance = AveExpr
      ) |>
      enrichVolcano::as_da(species = "Mus musculus") |>
      fgsea_of(gs$tested[[u]], TEST_WINDOW[[u]]) |>
      mutate(model = model, sex = sex, universe = u, .before = 1)
  }) |>
  list_rbind() |>
  mark_shown(c("model", "sex", "contrast")) |>
  transmute(
    model, sex, universe,
    collection = database, phenotype = contrast, term, NES = score, fgsea_p = p,
    fgsea_fdr = padj, size, overlap_group, shown, leading_edge = flatten(leading_edge)
  ) |>
  left_join(
    fit$score_associations |>
      filter(model %in% association_models, level == "pathway") |>
      select(model, sex, phenotype,
        term = feature, score_slope = estimate, score_r = r,
        score_r_veh = r_veh, score_r_pas = r_pas, score_p = P.Value, score_fdr = adj.P.Val
      ),
    by = c("model", "sex", "phenotype", "term")
  ) |>
  arrange(model, sex, phenotype, collection, fgsea_p)

stringdb <- string_db()
proteome <- filter(fit$protein_contrasts, universe == "full", !is.na(gene), gene != "")
string_map <- stringdb$map(
  data.frame(gene = unique(proteome$gene)), "gene",
  removeUnmappedRows = TRUE
)
stringdb$set_background(unique(string_map$STRING_id))
string_tests <- map(COMPARISON_CONTRASTS, \(cn) {
  hits <- filter(proteome, contrast == cn)
  lists <- list(
    `Π ≤ 0.05, up` = hits$gene[hits$pi_score <= ALPHA & hits$logFC > 0],
    `Π ≤ 0.05, down` = hits$gene[hits$pi_score <= ALPHA & hits$logFC < 0],
    `FDR ≤ 0.05` = hits$gene[hits$adj.P.Val <= ALPHA]
  )
  imap(lists, \(genes, name) {
    ids <- unique(string_map$STRING_id[string_map$gene %in% genes])
    e <- if (length(ids) >= 2) {
      stringdb$get_ppi_enrichment(ids)
    } else {
      list(edges = 0, lambda = 0, enrichment = 1)
    }
    tibble(
      contrast = cn, list = name, n_proteins = length(genes), mapped = length(ids), edges = e$edges,
      expected = e$lambda, fold = if_else(e$lambda > 0, e$edges / e$lambda, NA_real_),
      p = e$enrichment
    )
  }) |>
    list_rbind()
}) |>
  list_rbind() |>
  mutate(fdr = p.adjust(p, "BH"), .by = contrast)

# camera beside every fgsea result: the same competitive question tested across mice, with the
# set's own inter-protein correlation estimated (inter.gene.cor = NA) instead of assumed away (Wu &
# Smyth 2012). Unweighted: with a weight per protein, limma fits each protein's QR separately, the
# residual bases differ by rounding noise, and the estimated correlation, which compares residuals
# across proteins, moved male p-values by up to 45% under 1e-15 changes in the weights. Each set
# indexes the proteins of its member genes.
expression <- map(fit$E_norm, \(m) m[, fit$targets$sample_id])
set_index <- imap(sets, \(by_set, u) {
  gene_of <- deframe(distinct(filter(fit$protein_contrasts, universe == u), feature, gene))
  gene_of <- gene_of[rownames(expression[[u]])]
  map(by_set, \(members) which(gene_of %in% members))
})
camera_of <- function(u, design, contrast) {
  limma::camera(expression[[u]], set_index[[u]], design,
    contrast = contrast, inter.gene.cor = NA, sort = FALSE
  ) |>
    as_tibble(rownames = "term") |>
    transmute(universe = u, term, camera_p = PValue)
}
camera_results <- expand_grid(contrast = COMPARISON_CONTRASTS, u = names(gs$tested)) |>
  pmap(\(contrast, u) {
    camera_of(u, d$design, d$contrasts[, contrast]) |>
      mutate(contrast = contrast, .before = 1)
  }) |>
  list_rbind()
results <- results |>
  left_join(camera_results, by = c("contrast", "universe", "term")) |>
  mutate(
    camera_fdr = p.adjust(camera_p, "BH"),
    .by = c(contrast, universe, collection), .after = camera_p
  )
# For an association, camera tests the phenotype slope within one sex on that model's design.
slope_columns <- c(female = "slope_F", male = "slope_M")
association_camera <- fit$association_designs |>
  pmap(\(phenotype, model, design) {
    expand_grid(sex = names(slope_columns), u = names(gs$tested)) |>
      pmap(\(sex, u) {
        camera_of(u, design, which(colnames(design) == slope_columns[[sex]])) |>
          mutate(model = model, sex = sex, phenotype = phenotype, .before = 1)
      }) |>
      list_rbind()
  }) |>
  list_rbind()
association_fgsea <- association_fgsea |>
  left_join(association_camera, by = c("model", "sex", "phenotype", "universe", "term")) |>
  mutate(
    camera_fdr = p.adjust(camera_p, "BH"),
    .by = c(model, sex, phenotype, universe, collection)
  ) |>
  relocate(camera_p, camera_fdr, .after = fgsea_fdr)

# Does one sex's PAS response carry into the other's? The up and down proteins at Π ≤ 0.05 in one
# sex (full proteome) are tested in the other sex's PAS contrast by limma::fry, a self-contained
# test across mice, on the design with that contrast as its first coefficient. Unweighted and
# unstandardised, for the reason camera is unweighted; fry's standardisation also drops each
# protein's largest residual component, which depends on the order the mice are listed in. The
# share of members whose fold change in the other sex points the same way is reported beside it;
# members are listed by gene for the barcode, which places them in the other sex's ranking.
signature_proteins <- imap(SEX_CONTRASTS, \(cn, sex) {
  hits <- filter(proteome, contrast == cn, pi_score <= ALPHA)
  list(up = hits$feature[hits$logFC > 0], down = hits$feature[hits$logFC < 0])
})
signatures <- imap(SEX_CONTRASTS, \(cn, source) {
  target <- setdiff(names(SEX_CONTRASTS), source)
  lists <- signature_proteins[[source]]
  names(lists) <- paste(source, names(lists))
  design <- limma::contrastAsCoef(d$design, d$contrasts[, SEX_CONTRASTS[[target]]])$design
  tested <- limma::fry(expression$full, lists, design, contrast = 1, standardize = "none")
  in_target <- filter(proteome, contrast == SEX_CONTRASTS[[target]])
  imap(lists, \(proteins, name) {
    target_fc <- in_target$logFC[match(proteins, in_target$feature)]
    up <- str_ends(name, "up")
    members <- intersect(
      unique(in_target$gene[match(proteins, in_target$feature)]),
      names(ranks[[SEX_CONTRASTS[[target]]]]$full)
    )
    tibble(
      source = source, ranked_in = target, direction = str_remove(name, "^\\w+ "),
      n_proteins = length(proteins), n_ranked = length(members),
      same_way = mean(if (up) target_fc > 0 else target_fc < 0, na.rm = TRUE),
      fry_direction = tested[name, "Direction"], p = tested[name, "PValue"],
      p_mixed = tested[name, "PValue.Mixed"], members = paste(members, collapse = ";")
    )
  }) |>
    list_rbind()
}) |>
  list_rbind()

# What a list of proteins shares: ORA (fgsea::fora) of its genes against every detected gene, on
# the full-proteome collections; BH within collection. Descriptive: every list here comes from
# proteins chosen for responding.
full_sets <- gs$tested$full |>
  imap(\(sets, db) tibble(collection = db, term = names(sets), genes = sets)) |>
  list_rbind()
detected_genes <- names(ranks[[1]]$full)
ora_of <- function(genes) {
  map(split(full_sets, full_sets$collection), \(sets) {
    fgsea::fora(set_names(sets$genes, sets$term), genes, detected_genes,
      minSize = TEST_WINDOW$full[1], maxSize = TEST_WINDOW$full[2]
    ) |>
      as_tibble() |>
      transmute(
        collection = sets$collection[1], term = pathway, overlap, size, p = pval, fdr = padj,
        genes_in_list = length(genes),
        overlap_genes = map_chr(overlapGenes, \(g) paste(g, collapse = ";"))
      )
  }) |>
    list_rbind()
}

# The S3 heatmap rows: each comparison contrast's Π <= 0.05 proteins, split by direction, and each
# direction block cut into co-expression clusters by dynamicTreeCut's hybrid method on
# coexpression_tree() of the drawn values (Langfelder et al. 2008). Rows are numbered in drawing
# order: block by block, tree order within a cluster, clusters C1, C2, ... down the panel, and the
# proteins the cut leaves unassigned (label 0) at their block's end with no cluster.
CLUSTER_CUT <- list(deep_split = 2L, min_size = 10L)
targets <- arrange(fit$targets, factor(group, FIGURE_GROUPS))
observations <- readRDS(inputs[["proteins"]])$proteins$other$n.observations
block_clusters <- function(x) {
  tree <- coexpression_tree(x)
  cluster <- integer(nrow(x))
  cluster[tree$on_tree] <- dynamicTreeCut::cutreeDynamic(tree$tree,
    distM = tree$distance, method = "hybrid", deepSplit = CLUSTER_CUT$deep_split,
    minClusterSize = CLUSTER_CUT$min_size, verbose = 0
  )
  drawn <- tree$order[order(cluster[tree$order] == 0)]
  tibble(index = drawn, cut = cluster[drawn])
}
heatmap_clusters <- map(set_names(COMPARISON_CONTRASTS), \(cn) {
  selected <- fit$protein_contrasts |>
    filter(universe == "full", contrast == cn, pi_score <= ALPHA) |>
    mutate(direction = if_else(logFC > 0, "up", "down"))
  x <- heatmap_values(
    fit$E_norm$full[selected$feature, targets$sample_id], targets,
    observations[selected$feature, targets$sample_id] > 0
  )
  map(c("up", "down"), \(d) {
    rows <- which(selected$direction == d)
    mutate(block_clusters(x[rows, , drop = FALSE]), index = rows[index])
  }) |>
    list_rbind() |>
    mutate(
      block_cut = if_else(cut == 0, NA, paste(selected$direction[index], cut)),
      cluster = if_else(
        is.na(block_cut), NA, paste0("C", match(block_cut, unique(na.omit(block_cut))))
      )
    ) |>
    transmute(
      contrast = cn, row = row_number(), direction = selected$direction[index], cluster,
      feature = selected$feature[index], gene = selected$gene[index]
    )
}) |>
  list_rbind()
# ORA of each cluster; BH within cluster and collection.
cluster_ora <- heatmap_clusters |>
  filter(!is.na(cluster), !is.na(gene), gene != "") |>
  reframe(ora_of(unique(gene)), .by = c(contrast, direction, cluster)) |>
  rename(genes_in_cluster = genes_in_list) |>
  arrange(
    factor(contrast, COMPARISON_CONTRASTS), as.integer(str_remove(cluster, "C")), collection, p
  )

tier_table <- results |>
  reframe(tier_counts(p, padj), .by = c(contrast, universe, collection))
association_summary <- association_fgsea |>
  reframe(tier_counts(fgsea_p, fgsea_fdr),
    camera_fdr_le_0.05 = sum(camera_fdr <= ALPHA),
    singscore_p_le_0.05 = sum(score_p <= ALPHA),
    singscore_fdr_le_0.05 = sum(score_fdr <= ALPHA),
    .by = c(model, sex, universe, collection, phenotype)
  )

plot_index <- preview_pdf(
  file.path(out$reports, "02_enrichment.pdf"),
  set_names(
    map(COMPARISON_CONTRASTS, \(cn) {
      full <- filter(results, contrast == cn, universe == "full")
      \() gsea_table_page(full, sets$full, ranks[[cn]]$full)
    }),
    paste("strongest sets,", CONTRAST_LABELS[COMPARISON_CONTRASTS])
  )
)

saveRDS(
  list(
    results = results, ranks = ranks, sets = sets, compartments = compartments,
    literature = literature, association_fgsea = association_fgsea,
    signatures = signatures, heatmap_clusters = heatmap_clusters, cluster_ora = cluster_ora,
    string = string_tests, summary = tier_table
  ),
  file.path(out$data, "enrichment.rds")
)
write_workbook(
  list(
    summary = tier_table,
    fgsea = rename(results, set = term, fdr = padj),
    compartments = rename(compartments, fdr = padj),
    literature = rename(literature, fdr = list_fdr),
    association_summary = association_summary,
    association_fgsea = rename(association_fgsea, set = term, p = fgsea_p, fdr = fgsea_fdr),
    signatures = rename(signatures, size = n_proteins),
    cluster_ora = cluster_ora |>
      filter(p <= ALPHA) |>
      rename(set = term, members = overlap_genes),
    string = rename(string_tests, size = n_proteins), plot_index = plot_index
  ),
  c(
    "Per contrast and collection: sets tested, at p ≤ 0.05, FDR tiers.",
    paste(
      "Every tested set and contrast: NES, p, BH FDR, size, overlap group and whether F03 shows",
      "the set (one per group of overlapping sets at FDR ≤ 0.05), camera p and FDR",
      "(inter-protein correlation estimated), singscore contrast (log2 FC, p, FDR, r with",
      "group), GO Slim families, leading edge."
    ),
    str_glue(
      "GO Slim cellular compartments, fgsea per contrast, full proteome; compartments with more ",
      "than {SIDE_LIST_SIZE$compartments[2]} of the {length(ranks[[1]]$full)} ranked genes are ",
      "not tested."
    ),
    "Literature-informed and ROS sets, fgsea per contrast, full proteome, BH across the list.",
    "Per model, sex, collection and phenotype: fgsea, camera and singscore association counts.",
    paste(
      "Every set's fgsea on protein association t per model (pooled, treatment), sex and",
      "phenotype, with its overlap group and shown flag (picked as in fgsea, both universes",
      "together), camera p and FDR on the same slope, beside its singscore association (slope,",
      "r, p, FDR)."
    ),
    paste(
      "Each sex's Π ≤ 0.05 up and down proteins tested in the other sex's PAS contrast by fry:",
      "size, share moving the same way, fry direction, p and mixed p."
    ),
    paste(
      "ORA of each S3 heatmap co-expression cluster, full-proteome collections, sets at p ≤ 0.05:",
      "overlap, size, p, BH FDR within cluster and collection (over every set tested), members",
      "of the set in the cluster. The cluster of each protein is in enrichment.rds",
      "(heatmap_clusters)."
    ),
    "STRING v12 connectivity of each contrast's Π and FDR lists.",
    "Which PDF page shows what."
  ),
  out$workbook
)

tier_table
