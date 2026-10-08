# Every pathway test the comparisons read, run once. fgsea (enrichVolcano::run_enrichment) on the
# moderated t of each comparison's contrast: GO:BP and Reactome in the full proteome, MitoCarta in
# the mitochondrial universe, BH within contrast and collection over every tested set, with the same
# set's singscore contrast from 01_Fit beside it. GO Slim compartments, and the literature-informed
# and ROS sets as one list in the full proteome. Every term carries an EnrichmentMap redundancy flag
# (enrichVolcano::dedup_terms, similarity >= 0.375, the stronger term kept) that only picks display
# terms; nothing is dropped. fgsea also runs on each protein-phenotype association t per sex, beside
# the same set's singscore association. Cross-sex signatures: each sex's proteins at p <= 0.05, up
# and down, tested by fgsea in the other sex's ranking. STRING v12 connectivity of each contrast's Π
# and FDR lists is descriptive: subunits of one complex vary together, so connected lists are
# expected.

source(here::here("R", "helpers.R"))

out <- stage_paths("02_Statistics", "02_Enrichment")
study_files <- expand_grid(
  universe = c("full", "mito"), file = c("da_results", "matrix", "weights", "samples", "contrasts")
) |>
  with(set_names(
    here("02_Statistics", "01_Fit", "c_data", paste0("study_", universe), paste0(file, ".csv")),
    paste0("study_", universe, "_", file)
  ))
inputs <- c(
  gene_sets = here("01_Processing", "04_Gene_Sets", "c_data", "gene_sets.rds"),
  fit = here("02_Statistics", "01_Fit", "c_data", "fit.rds"),
  study_files
)
gs <- readRDS(inputs[["gene_sets"]])
fit <- readRDS(inputs[["fit"]])
SIMILARITY_CUTOFF <- 0.375

sets <- map(gs$tested, \(x) unlist(unname(x), recursive = FALSE))
studies <- map(set_names(names(gs$tested)), \(u) {
  suppressMessages(enrichVolcano::read_study(dirname(inputs[[paste0("study_", u, "_da_results")]]),
    species = "Mus musculus"
  ))
})

# One contrast at a time and seeded each time, so a contrast's result does not depend on which
# others ran beside it.
fgsea_of <- function(study, collection, contrast, size) {
  study$da <- study$da[study$da$contrast == contrast, ]
  study$contrasts <- study$contrasts[study$contrasts$name == contrast, , drop = FALSE]
  set.seed(SEED)
  suppressMessages(enrichVolcano::run_enrichment(study, collection,
    tests = "fgsea", min_size = size[1], max_size = size[2]
  ))$fgsea
}
# The ranking fgsea tested: moderated t, one protein per gene by abundance.
gene_ranks <- function(da) {
  keep <- one_per_gene(da$protein, da$gene, da$abundance)
  set_names(da$rank[keep], da$gene[keep])
}
flatten <- \(x) map_chr(x, \(g) paste(g, collapse = ";"))

enrichment_objects <- map(set_names(COMPARISON_CONTRASTS), \(cn) {
  imap(gs$tested, \(collection, u) {
    enrichVolcano::dedup_terms(fgsea_of(studies[[u]], collection, cn, TEST_WINDOW[[u]]), sets[[u]],
      cutoff = SIMILARITY_CUTOFF, term_threshold = Inf
    )
  })
})
ranks <- map(set_names(COMPARISON_CONTRASTS), \(cn) {
  map(studies, \(s) gene_ranks(s$da[s$da$contrast == cn, ]))
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
stopifnot(!anyNA(results$score_p), !anyNA(results$dedup_status))

side_list <- function(collection, size) {
  map(COMPARISON_CONTRASTS, \(cn) fgsea_of(studies$full, collection, cn, size)@results) |>
    list_rbind() |>
    transmute(
      contrast,
      set = term, NES = score, p, padj, size, direction, leading_edge = flatten(leading_edge)
    )
}
compartments <- side_list(gs$compartments, c(15, Inf))
literature <- side_list(list(Literature = gs$literature_sets), c(5, Inf)) |>
  rename(list_fdr = padj) |>
  left_join(gs$literature, by = "set") |>
  relocate(hypothesis, expected, .after = contrast)

# Pathway associations: fgsea on each protein's association t with one phenotype across one sex's
# PAS and vehicle mice (genes permuted, so its FDR runs optimistic), beside the same set's singscore
# association, which is tested across mice.
association_fgsea <- expand_grid(sex = c("female", "male"), u = names(gs$tested)) |>
  pmap(\(sex, u) {
    da <- fit$protein_associations |>
      filter(!adjusted, sex == !!sex, universe == u) |>
      transmute(
        protein = feature, gene, contrast = phenotype, logFC = estimate, t, P.Value, adj.P.Val,
        abundance = AveExpr
      ) |>
      enrichVolcano::as_da(species = "Mus musculus")
    set.seed(SEED)
    suppressMessages(enrichVolcano::run_enrichment(da, gs$tested[[u]],
      tests = "fgsea", min_size = TEST_WINDOW[[u]][1], max_size = TEST_WINDOW[[u]][2]
    ))$fgsea@results |>
      mutate(sex = sex, universe = u, .before = 1)
  }) |>
  list_rbind() |>
  transmute(
    sex, universe,
    collection = database, phenotype = contrast, term, NES = score, fgsea_p = p,
    fgsea_fdr = padj, size, leading_edge = flatten(leading_edge)
  ) |>
  left_join(
    fit$score_associations |>
      filter(!adjusted, level == "pathway") |>
      select(sex, phenotype,
        term = feature, score_slope = estimate, score_r = r,
        score_r_veh = r_veh, score_r_pas = r_pas, score_p = P.Value, score_fdr = adj.P.Val
      ),
    by = c("sex", "phenotype", "term")
  ) |>
  arrange(sex, phenotype, collection, fgsea_p)
stopifnot(!anyNA(association_fgsea$score_p))

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

# Does one sex's PAS response carry into the other's? The up and down proteins at p <= 0.05 in one
# sex (one per gene, full proteome) are placed in the other sex's ranking; fgsea tests each list
# there, and the share of members ranked on the same side as in the source sex is reported.
SEXES <- c(female = "PAS_vs_VEH_female", male = "PAS_vs_VEH_male")
signature_genes <- imap(SEXES, \(cn, sex) {
  hits <- filter(proteome, contrast == cn, P.Value <= ALPHA)
  list(up = unique(hits$gene[hits$logFC > 0]), down = unique(hits$gene[hits$logFC < 0]))
})
signatures <- imap(SEXES, \(cn, source) {
  target <- setdiff(names(SEXES), source)
  in_target <- ranks[[SEXES[[target]]]]$full
  lists <- set_names(signature_genes[[source]], paste(source, names(signature_genes[[source]])))
  tested <- fgsea_of(studies$full, list(Signature = lists), SEXES[[target]], c(5, Inf))@results
  imap(lists, \(genes, name) {
    ranked <- in_target[intersect(genes, names(in_target))]
    row <- tested[tested$term == name, ]
    tibble(
      source = source, ranked_in = target, direction = str_remove(name, "^\\w+ "),
      n_proteins = length(genes), n_ranked = length(ranked),
      same_way = mean(if (str_ends(name, "up")) ranked > 0 else ranked < 0),
      NES = row$score, p = row$p, members = paste(names(ranked), collapse = ";")
    )
  }) |>
    list_rbind()
}) |>
  list_rbind()
stopifnot(nrow(signatures) == 4, !anyNA(signatures$p))

tier_table <- results |>
  reframe(tier_counts(p, padj), .by = c(contrast, universe, collection)) |>
  select(-pi_le_0.05)
association_summary <- association_fgsea |>
  reframe(tier_counts(fgsea_p, fgsea_fdr),
    singscore_p_le_0.05 = sum(score_p <= ALPHA),
    singscore_fdr_le_0.05 = sum(score_fdr <= ALPHA), .by = c(sex, universe, collection, phenotype)
  ) |>
  select(-pi_le_0.05)

summary_page <- tier_table |>
  select(contrast, collection, tests, p_le_0.05, expected_by_chance, fdr_le_0.05) |>
  pivot_longer(c(p_le_0.05, fdr_le_0.05), names_to = "tier", values_to = "n") |>
  mutate(
    contrast = factor(CONTRAST_LABELS[contrast], CONTRAST_LABELS[COMPARISON_CONTRASTS]),
    tier = factor(tier, c("p_le_0.05", "fdr_le_0.05"), c("p ≤ 0.05", "FDR ≤ 0.05")),
    collection = factor(collection, c("GO:BP", "Reactome", "MitoCarta"))
  ) |>
  ggplot(aes(n, collection, fill = tier)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_point(
    aes(x = expected_by_chance),
    data = \(d) filter(d, tier == "p ≤ 0.05"), shape = 124, size = 3,
    position = position_nudge(y = -0.2)
  ) +
  geom_text(aes(label = n), position = position_dodge(width = 0.8), hjust = -0.2, size = 2.2) +
  facet_wrap(~contrast, ncol = 1) +
  scale_fill_manual(values = c("p ≤ 0.05" = "grey65", "FDR ≤ 0.05" = "black"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "Sets", y = NULL, title = "fgsea, each comparison's contrast",
    subtitle = "BH within contrast and collection; bar mark: sets at p ≤ 0.05 expected by chance"
  ) +
  theme_figure()
association_page <- association_summary |>
  mutate(
    sex = if_else(sex == "female", "Female", "Male"),
    trait = factor(PHENOTYPE_LABELS[phenotype], rev(PHENOTYPE_LABELS)),
    collection = factor(collection, c("GO:BP", "Reactome", "MitoCarta"))
  ) |>
  ggplot(aes(collection, trait, fill = p_le_0.05 / expected_by_chance)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = str_glue("{p_le_0.05} / {tests}\nFDR {fdr_le_0.05}")), size = 1.9) +
  facet_wrap(~sex) +
  scale_fill_gradient(low = "white", high = "grey35", name = "p ≤ 0.05 /\nchance") +
  labs(
    x = NULL, y = NULL, title = "Pathways moving with each phenotype",
    subtitle = "fgsea on protein association t across one sex's PAS and vehicle mice (co-movement)"
  ) +
  theme_figure()
string_page <- string_tests |>
  filter(expected > 0) |>
  mutate(contrast = factor(CONTRAST_LABELS[contrast], CONTRAST_LABELS[COMPARISON_CONTRASTS])) |>
  ggplot(aes(fold, list)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(size = n_proteins), colour = "grey30") +
  geom_text(
    aes(label = sprintf("%d edges, p %s", edges, if_else(p == 0, "< 1e-16", sprintf("%.2g", p)))),
    nudge_y = 0.3, size = 2
  ) +
  facet_wrap(~contrast, ncol = 1) +
  scale_size_area(max_size = 4, name = "Proteins") +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  labs(
    x = "STRING edges / expected", y = NULL, title = "STRING connectivity of each protein list",
    subtitle = paste(
      "v12, score ≥ 400, against the detected proteome; descriptive;",
      "lists with no expected edges omitted"
    )
  ) +
  theme_figure()
plot_index <- write_report(file.path(out$reports, "02_enrichment.pdf"), list(
  "fgsea sets per tier, each contrast" = summary_page,
  "pathway associations per phenotype" = association_page,
  "STRING connectivity" = string_page
))

saveRDS(
  list(
    results = results, objects = enrichment_objects, ranks = ranks, da = map(studies, "da"),
    sets = sets,
    compartments = compartments, literature = literature, association_fgsea = association_fgsea,
    signatures = signatures, string = string_tests, summary = tier_table
  ),
  file.path(out$data, "enrichment.rds")
)
write_workbook(
  list(
    summary = tier_table, fgsea = results, compartments = compartments, literature = literature,
    association_summary = association_summary, association_fgsea = association_fgsea,
    signatures = signatures, string = string_tests, plot_index = plot_index
  ),
  c(
    "Per contrast and collection: sets tested, at p ≤ 0.05, expected by chance, FDR tiers.",
    paste(
      "Every tested set and contrast: NES, p, BH FDR, size, redundancy flag, singscore",
      "contrast (log2 FC, p, FDR, r with group), GO Slim families, leading edge."
    ),
    "GO Slim cellular compartments, fgsea per contrast, full proteome.",
    "Literature-informed and ROS sets, fgsea per contrast, full proteome, BH across the list.",
    "Per sex, collection and phenotype: fgsea and singscore association counts.",
    paste(
      "Every set's fgsea on protein association t per sex and phenotype,",
      "beside its singscore association (slope, r, p, FDR)."
    ),
    paste(
      "Each sex's p ≤ 0.05 up and down proteins in the other sex's ranking:",
      "size, share ranked the same way, fgsea NES and p."
    ),
    "STRING v12 connectivity of each contrast's Π and FDR lists.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

tier_table

sessionInfo()
