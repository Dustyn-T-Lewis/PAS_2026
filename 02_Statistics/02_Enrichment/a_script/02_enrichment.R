# Every pathway test the comparisons read, run once. fgsea (enrichVolcano::run_enrichment) on the
# moderated t of each comparison's contrast: GO:BP, Reactome, Hallmark and KEGG in the full
# proteome, MitoCarta in the mitochondrial universe, BH within contrast and collection over every
# tested set, with the same set's singscore contrast from 01_Fit beside it. GO Slim compartments,
# and the literature-informed and ROS sets as one list in the full proteome. Every FDR-level term
# carries a redundancy flag from fgsea::collapsePathways (through enrichVolcano::dedup_terms):
# a term is redundant when it is no longer enriched once a stronger term's genes are conditioned
# on. The flag only picks display terms; nothing is dropped. fgsea also runs on each
# protein-phenotype association t per sex and association model. camera sits beside every fgsea
# result, contrast and association, testing the same sets across mice. Cross-sex signatures: each
# sex's proteins at Π ≤ 0.05, up and down, tested by limma::fry in the other sex's contrast
# across mice. ORA (fgsea::fora) describes each response pattern of the PAS F-test proteins.
# STRING v12 connectivity of each contrast's Π and FDR lists is descriptive: subunits of one
# complex vary together.

source(here::here("00_Setup", "Shared_functions", "helpers.R"))

out <- stage_paths("02_Statistics", "02_Enrichment")
inputs <- c(
  gene_sets = output_path("01_Processing", "04_Gene_Sets", "gene_sets.rds"),
  design = output_path("01_Processing", "03_Design", "design.rds"),
  fit = output_path("02_Statistics", "01_Fit", "fit.rds"),
  STRING_FILES
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

# One contrast at a time and seeded each time, so a contrast's result does not depend on which
# others ran beside it.
fgsea_of <- function(da, collection, size) {
  set.seed(SEED)
  suppressMessages(enrichVolcano::run_enrichment(da, collection,
    tests = "fgsea", min_size = size[1], max_size = size[2]
  ))$fgsea
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
enrichment_objects <- map(set_names(COMPARISON_CONTRASTS), \(cn) {
  imap(gs$tested, \(collection, u) {
    tested <- fgsea_of(da[[u]][da[[u]]$contrast == cn, ], collection, TEST_WINDOW[[u]])
    set.seed(SEED)
    enrichVolcano::dedup_terms(tested, sets[[u]],
      method = "collapse_pathways", term_threshold = ALPHA,
      stats = set_names(list(ranks[[cn]][[u]]), cn)
    )
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
results <- imap(enrichment_objects, \(by_universe, cn) {
  imap(by_universe, \(e, u) mutate(e@results, universe = u, .before = 1)) |>
    list_rbind()
}) |>
  list_rbind() |>
  left_join(
    summarise(families, families = paste(family, collapse = "; "), .by = term),
    by = "term"
  ) |>
  left_join(score_stats, by = c("contrast", "term")) |>
  transmute(
    contrast, universe,
    collection = database, term, NES = score, p, padj, size, direction,
    dedup_status, merged_into, similarity, score_logFC, score_p, score_fdr, r_group, families,
    leading_edge = flatten(leading_edge)
  ) |>
  arrange(factor(contrast, COMPARISON_CONTRASTS), universe, collection, p)

# The side lists run on fgsea directly, on the ranking the main tests used: enrichVolcano keeps
# fgsea's default nPermSimple, and at that default a set with unbalanced gene-level statistics
# gets no p-value and drops out. fgsea's own advice is nPermSimple = 10000.
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
compartments <- side_list(gs$compartments, c(15, Inf))
literature <- side_list(list(Literature = gs$literature_sets), c(5, Inf)) |>
  rename(list_fdr = padj) |>
  left_join(gs$literature, by = "set") |>
  relocate(hypothesis, expected, .after = contrast)

# Pathway associations: fgsea on each protein's association t with one phenotype within one sex,
# under both association models (genes permuted, so its FDR runs optimistic), beside the same set's
# singscore association under the same model, which is tested across mice. Each phenotype column
# carries the same collapsePathways flag as the contrasts.
association_models <- unique(fit$association_designs$model)
association_fgsea <- expand_grid(
  model = association_models, sex = names(SEX_COMPARISONS), u = names(gs$tested)
) |>
  pmap(\(model, sex, u) {
    da <- fit$protein_associations |>
      filter(model == !!model, sex == !!sex, universe == u) |>
      transmute(
        protein = feature, gene, contrast = phenotype, logFC = estimate, t, P.Value, adj.P.Val,
        abundance = AveExpr
      ) |>
      enrichVolcano::as_da(species = "Mus musculus")
    set.seed(SEED)
    tested <- suppressMessages(enrichVolcano::run_enrichment(da, gs$tested[[u]],
      tests = "fgsea", min_size = TEST_WINDOW[[u]][1], max_size = TEST_WINDOW[[u]][2]
    ))$fgsea
    set.seed(SEED)
    enrichVolcano::dedup_terms(tested, sets[[u]],
      method = "collapse_pathways", term_threshold = ALPHA,
      stats = map(split(da, da$contrast), gene_ranks)
    )@results |>
      mutate(model = model, sex = sex, universe = u, .before = 1)
  }) |>
  list_rbind() |>
  transmute(
    model, sex, universe,
    collection = database, phenotype = contrast, term, NES = score, fgsea_p = p,
    fgsea_fdr = padj, size, dedup_status, merged_into, leading_edge = flatten(leading_edge)
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

# Does one sex's PAS response carry into the other's? The up and down proteins at Π ≤ 0.05 in one
# sex (full proteome) are tested in the other sex's PAS contrast by limma::fry, a self-contained
# test across mice, on the design with that contrast as its first coefficient. Unweighted and
# unstandardised, for the reason camera is unweighted; fry's standardisation also drops each
# protein's largest residual component, which depends on the order the mice are listed in. The
# share of members whose fold change in the other sex points the same way is reported beside it;
# members are listed by gene for the barcode, which places them in the other sex's ranking.
SEXES <- map_chr(SEX_COMPARISONS, "contrast")
full_matrix <- fit$E_norm$full[, fit$targets$sample_id]

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

signature_proteins <- imap(SEXES, \(cn, sex) {
  hits <- filter(proteome, contrast == cn, pi_score <= ALPHA)
  list(up = hits$feature[hits$logFC > 0], down = hits$feature[hits$logFC < 0])
})
signatures <- imap(SEXES, \(cn, source) {
  target <- setdiff(names(SEXES), source)
  lists <- signature_proteins[[source]]
  names(lists) <- paste(source, names(lists))
  design <- limma::contrastAsCoef(d$design, d$contrasts[, SEXES[[target]]])$design
  tested <- limma::fry(full_matrix, lists, design, contrast = 1, standardize = "none")
  in_target <- filter(proteome, contrast == SEXES[[target]])
  imap(lists, \(proteins, name) {
    target_fc <- in_target$logFC[match(proteins, in_target$feature)]
    up <- str_ends(name, "up")
    members <- intersect(
      unique(in_target$gene[match(proteins, in_target$feature)]),
      names(ranks[[SEXES[[target]]]]$full)
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

# What the PAS F-test proteins share: ORA (fgsea::fora) of each response pattern's genes against
# every detected gene, on the full-proteome collections; BH within pattern and collection.
# Descriptive: the patterns come from proteins chosen for responding.
full_sets <- gs$tested$full |>
  imap(\(sets, db) tibble(collection = db, term = names(sets), genes = sets)) |>
  list_rbind()
detected_genes <- names(ranks[[1]]$full)
pattern_genes <- fit$protein_ftest |>
  filter(universe == "full", P.Value <= ALPHA, !is.na(gene), gene != "") |>
  summarise(genes = list(unique(gene)), proteins = n(), .by = pattern)
pattern_ora <- pmap(pattern_genes, \(pattern, genes, proteins) {
  map(split(full_sets, full_sets$collection), \(sets) {
    fgsea::fora(set_names(sets$genes, sets$term), genes, detected_genes,
      minSize = TEST_WINDOW$full[1], maxSize = TEST_WINDOW$full[2]
    ) |>
      as_tibble() |>
      transmute(
        pattern = pattern, collection = sets$collection[1], term = pathway, overlap, size,
        p = pval, fdr = p.adjust(pval, "BH"), genes_in_pattern = length(genes),
        overlap_genes = map_chr(overlapGenes, \(g) paste(g, collapse = ";"))
      )
  }) |>
    list_rbind()
}) |>
  list_rbind() |>
  arrange(pattern, collection, p)

tier_table <- results |>
  reframe(tier_counts(p, padj), .by = c(contrast, universe, collection)) |>
  select(-pi_le_0.05)
association_summary <- association_fgsea |>
  reframe(tier_counts(fgsea_p, fgsea_fdr),
    camera_fdr_le_0.05 = sum(camera_fdr <= ALPHA),
    singscore_p_le_0.05 = sum(score_p <= ALPHA),
    singscore_fdr_le_0.05 = sum(score_fdr <= ALPHA),
    .by = c(model, sex, universe, collection, phenotype)
  ) |>
  select(-pi_le_0.05)

# fgsea's own table of running scores: the five strongest sets each way per comparison contrast.
gsea_table <- \(cn) {
  \() {
    top <- results |>
      filter(contrast == cn, universe == "full") |>
      slice_min(p, n = 5, by = direction, with_ties = FALSE) |>
      arrange(desc(NES)) |>
      transmute(pathway = term, NES, pval = p, padj)
    short <- str_trunc(enrichVolcano::clean_label(top$pathway, width = 200), 36)
    pathways <- set_names(sets$full[top$pathway], short)
    print(fgsea::plotGseaTable(pathways, ranks[[cn]]$full, mutate(top, pathway = short),
      gseaParam = 1
    ))
  }
}
plot_index <- preview_pdf(
  file.path(out$reports, "02_enrichment.pdf"),
  set_names(
    map(COMPARISON_CONTRASTS, gsea_table),
    paste("strongest sets,", CONTRAST_LABELS[COMPARISON_CONTRASTS])
  )
)

saveRDS(
  list(
    results = results, ranks = ranks, sets = sets, compartments = compartments,
    literature = literature, association_fgsea = association_fgsea,
    signatures = signatures, pattern_ora = pattern_ora, string = string_tests, summary = tier_table
  ),
  file.path(out$data, "enrichment.rds")
)
write_workbook(
  list(
    summary = tier_table, fgsea = results, compartments = compartments, literature = literature,
    association_summary = association_summary, association_fgsea = association_fgsea,
    signatures = signatures, pattern_ora = pattern_ora, string = string_tests,
    plot_index = plot_index
  ),
  c(
    "Per contrast and collection: sets tested, at p ≤ 0.05, expected by chance, FDR tiers.",
    paste(
      "Every tested set and contrast: NES, p, BH FDR, size, collapsePathways flag, camera p and",
      "FDR (inter-protein correlation estimated), singscore",
      "contrast (log2 FC, p, FDR, r with group), GO Slim families, leading edge."
    ),
    "GO Slim cellular compartments, fgsea per contrast, full proteome.",
    "Literature-informed and ROS sets, fgsea per contrast, full proteome, BH across the list.",
    "Per model, sex, collection and phenotype: fgsea, camera and singscore association counts.",
    paste(
      "Every set's fgsea on protein association t per model (pooled, treatment), sex and",
      "phenotype, with its collapsePathways flag and camera p and FDR on the same slope, beside",
      "its singscore association (slope, r, p, FDR)."
    ),
    paste(
      "Each sex's Π ≤ 0.05 up and down proteins tested in the other sex's PAS contrast by fry:",
      "size, share moving the same way, fry direction, p and mixed p."
    ),
    paste(
      "ORA of each response pattern of the PAS F-test proteins (p ≤ 0.05), full-proteome",
      "collections: overlap, size, p, BH FDR within pattern and collection, overlap genes."
    ),
    "STRING v12 connectivity of each contrast's Π and FDR lists.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

tier_table
