#!/usr/bin/env Rscript
# One singscore per mouse per tested GO:BP and Reactome set, CORUM complex and STRING cluster on
# the full matrix, and per tested MitoCarta pathway on the mitochondrial matrix. The GO:BP,
# Reactome and MitoCarta sets are the ones fgsea tests in 02_Enrichment. A score is the mean
# within-mouse rank of the group's members, so it uses no group label; 01_Fit tests it.

source(here::here("R", "helpers.R"))
library(singscore)

out <- stage_paths("01_Processing", "05_Scores")
inputs <- c(
  proteins = here("01_Processing", "02_Quantification", "c_data", "proteins.rds"),
  gene_sets = here("01_Processing", "04_Gene_Sets", "c_data", "gene_sets.rds")
)
quant <- readRDS(inputs[["proteins"]])
gs <- readRDS(inputs[["gene_sets"]])
targets <- quant$proteins$targets

# One protein per gene by average abundance, the rule enrichVolcano applies before fgsea.
symbol <- set_names(first_symbol(quant$proteins$genes$Genes), rownames(quant$proteins))
by_gene <- function(E) {
  keep <- one_per_gene(rownames(E), symbol[rownames(E)], rowMeans(E))
  m <- E[keep, ]
  rownames(m) <- symbol[rownames(m)]
  m
}
score <- function(E, collections) {
  sets <- unlist(unname(collections), recursive = FALSE)
  multiScore(rankGenes(by_gene(E)), upSetColc = sets)$Scores
}
scores <- list(
  full = score(quant$E_norm$full, c(gs$tested$full, gs$complexes, gs$clusters)),
  mito = score(quant$E_norm$mito, gs$tested$mito)
)
stopifnot(map_lgl(scores, \(s) identical(colnames(s), targets$sample_id)))

collection_of <- c(gs$tested$full, gs$tested$mito, gs$complexes, gs$clusters) |>
  imap(\(sets, db) set_names(rep(db, length(sets)), names(sets))) |>
  unname() |>
  unlist()
catalog <- imap(scores, \(s, universe) {
  tibble(
    universe = universe, set = rownames(s), collection = collection_of[rownames(s)],
    mean_score = rowMeans(s)
  )
}) |>
  list_rbind()

pages <- distinct(catalog, universe, collection) |>
  pmap(\(universe, collection) {
    rows <- catalog$set[catalog$universe == universe & catalog$collection == collection]
    group_heatmap(scores[[universe]][rows, , drop = FALSE], targets,
      title = str_glue("{collection} singscores"),
      subtitle = "Rows z-scored and clustered on 1 - Pearson r (display only); mice by group"
    )
  }) |>
  set_names(str_glue("singscore heatmap, {distinct(catalog, universe, collection)$collection}"))
plot_index <- write_report(file.path(out$reports, "05_scores.pdf"), pages)

saveRDS(list(scores = scores, collection = collection_of), file.path(out$data, "scores.rds"))
write_workbook(
  list(
    catalog = catalog,
    scores_full = as_tibble(scores$full, rownames = "set"),
    scores_mito = as_tibble(scores$mito, rownames = "set"),
    plot_index = plot_index
  ),
  c(
    "Every scored set, complex and cluster: universe, collection, mean score.",
    "Scores by mouse: GO:BP and Reactome sets, CORUM complexes and STRING clusters.",
    "Scores by mouse: MitoCarta pathways in the mitochondrial universe.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

count(catalog, universe, collection)

sessionInfo()
