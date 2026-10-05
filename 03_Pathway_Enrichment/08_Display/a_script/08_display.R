#!/usr/bin/env Rscript
# What the pathway figures draw, chosen once. The display pool (display_pool() in R/helpers.R) is
# collapsed with enrichVolcano's EnrichmentMap rule over every term and across databases; the 12
# strongest collapsed terms per sex, and for the treatment-by-sex contrast, per universe go on the
# rings. The term tree takes the 15
# strongest uncollapsed terms per sex in the full proteome, since grouping overlapping terms is its
# job. It uses the collapse's own similarity, the EnrichmentMap coefficient (mean of Jaccard and
# overlap coefficient; Merico et al. 2010), with average linkage, and cuts at the collapse's cut-off
# of 0.375, so a cluster holds terms as alike as those the rings merge; a term like no other stands
# alone. Each cluster is named after its member with the smallest p. The 10 strongest collapsed
# terms for each baseline sex contrast (male minus female, in vehicle and in PAS) go to Figure 1.
# The literature-informed sets (00_Gene_Sets) are taken in full with BH across the list within each
# contrast, so their FDR answers only for the sets named. Display only.

source(here::here("R", "helpers.R"))

SEXES <- c("PAS_vs_VEH_female", "PAS_vs_VEH_male")
SEX_DIFFERENCES <- c("sex_in_VEH", "sex_in_PAS")
N_SEX <- 10
N_RING <- 12
N_TREE <- 15
SIMILARITY_CUTOFF <- 0.375

out <- stage_paths("03_Pathway_Enrichment", "08_Display")
inputs <- c(
  gene_sets = here("03_Pathway_Enrichment", "00_Gene_Sets", "c_data", "gene_sets.rds"),
  enrichment = here("03_Pathway_Enrichment", "02_Contrasts", "c_data", "enrichment.rds")
)
gs <- readRDS(inputs[["gene_sets"]])
enrichment <- readRDS(inputs[["enrichment"]])

# dedup_terms() compares terms within one database, so the databases are pooled under one name for
# the call and restored after.
display <- imap(enrichment, \(tests, u) {
  e <- tests$fgsea
  e@results <- display_pool(e@results, gs$slim_terms)
  database <- e@results$database
  e@results$database <- "pooled"
  e <- enrichVolcano::dedup_terms(e, unlist(unname(gs[[u]]), recursive = FALSE),
    method = "enrichmentmap", term_threshold = 1
  )
  e@results$database <- database
  e
})

results <- list_rbind(imap(display, \(e, u) mutate(e@results, universe = u, .before = 1)))
pool <- filter(results, contrast %in% SEXES)
ring_terms <- results |>
  filter(contrast %in% c(SEXES, "treatment_by_sex"), dedup_status == "kept") |>
  slice_min(p, n = N_RING, by = c(universe, contrast), with_ties = FALSE)

tree_terms <- pool |>
  filter(universe == "full") |>
  slice_min(p, n = N_TREE, by = contrast, with_ties = FALSE) |>
  distinct(term) |>
  pull(term)
sets <- unlist(unname(gs$full), recursive = FALSE)[tree_terms]
similarity <- outer(tree_terms, tree_terms, Vectorize(\(a, b) {
  shared <- length(intersect(sets[[a]], sets[[b]]))
  mean(c(
    shared / length(union(sets[[a]], sets[[b]])),
    shared / min(length(sets[[a]]), length(sets[[b]]))
  ))
}))
dimnames(similarity) <- list(tree_terms, tree_terms)
tree <- hclust(as.dist(1 - similarity), method = "average")

clusters <- pool |>
  filter(universe == "full", term %in% tree_terms) |>
  summarise(p = min(p), .by = term) |>
  mutate(cluster = cutree(tree, h = 1 - SIMILARITY_CUTOFF)[term]) |>
  mutate(cluster_name = term[which.min(p)], size = n(), .by = cluster) |>
  select(term, cluster, cluster_name, size)
cohesion <- clusters |>
  filter(size > 1) |>
  summarise(
    terms = n(),
    within = mean(similarity[term, term][upper.tri(diag(n()))]),
    nearest_other = max(similarity[term, setdiff(tree_terms, term)]),
    .by = c(cluster, cluster_name)
  )
sex_terms <- results |>
  filter(universe == "full", contrast %in% SEX_DIFFERENCES, dedup_status == "kept") |>
  slice_min(p, n = N_SEX, by = contrast, with_ties = FALSE) |>
  distinct(term) |>
  pull(term)
sex_dots <- results |>
  filter(universe == "full", term %in% sex_terms, contrast %in% SEX_DIFFERENCES) |>
  select(contrast, database, term, NES = score, p, padj)
# Taken from every tested term, not the display pool: a named set is reported whatever it shows.
literature <- imap(enrichment, \(tests, u) mutate(tests$fgsea@results, universe = u)) |>
  list_rbind() |>
  filter(contrast %in% PRIMARY) |>
  inner_join(gs$literature, by = c(term = "set")) |>
  mutate(list_fdr = p.adjust(p, "BH"), .by = contrast) |>
  select(contrast, hypothesis, expected, database, term, NES = score, p, padj, list_fdr)
tree_dots <- display$full@results |>
  filter(term %in% tree_terms, contrast %in% CONTRASTS) |>
  select(contrast, database, term, NES = score, p, padj)

saveRDS(
  list(
    display = display, ring_terms = ring_terms, tree = tree, clusters = clusters,
    tree_dots = tree_dots, sex_dots = sex_dots, literature = literature
  ),
  file.path(out$data, "display.rds")
)
write_workbook(
  list(
    ring_terms = select(ring_terms, universe, contrast, database, term, NES = score, p, padj),
    tree_terms = left_join(clusters, tibble(term = tree$labels[tree$order], order = seq_along(
      tree$order
    )), by = "term") |> arrange(order),
    tree_dots = tree_dots, cohesion = cohesion, sex_dots = sex_dots,
    literature = literature
  ),
  c(
    "Ring terms: the 12 strongest collapsed terms per universe, for each sex and the interaction.",
    "Tree terms in tree order, with their cluster and the cluster's representative (smallest p).",
    "NES, p and FDR (fgsea, 02_Contrasts) for every tree term in every contrast.",
    "Per cluster: mean EnrichmentMap similarity inside it and its highest to any outside term.",
    "Baseline sex terms (male minus female; 10 strongest collapsed per contrast), both contrasts.",
    "Literature-informed sets per primary contrast: fgsea NES, p, FDR, and BH across the list."
  ),
  out$workbook, inputs
)

count(clusters, cluster, cluster_name)
cohesion
sessionInfo()
