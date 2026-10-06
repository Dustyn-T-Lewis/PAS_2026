#!/usr/bin/env Rscript
# One singscore per mouse per GO:BP term, CORUM complex and STRING cluster on the full matrix, and
# per MitoCarta pathway on the mitochondrial matrix. A score is the mean within-mouse rank of the
# group's members, so it uses no group label; 07_Fit tests it.

source(here::here("R", "helpers.R"))
library(singscore)

out <- stage_paths("00_Trunk", "05_Scores")
inputs <- c(
  proteins = here("00_Trunk", "02_Quantification", "c_data", "proteins.rds"),
  gene_sets = here("00_Trunk", "04_Gene_Sets", "c_data", "gene_sets.rds")
)
quant <- readRDS(inputs[["proteins"]])
gs <- readRDS(inputs[["gene_sets"]])
targets <- quant$proteins$targets

# One protein per gene by average abundance, the rule enrichVolcano applies before fgsea.
symbol <- set_names(str_split_i(quant$proteins$genes$Genes, ";", 1), rownames(quant$proteins))
by_gene <- function(E) {
  keep <- one_per_gene(rownames(E), symbol[rownames(E)], rowMeans(E))
  m <- E[keep, ]
  rownames(m) <- symbol[rownames(m)]
  m
}
score <- function(E, collections, min_size) {
  m <- by_gene(E)
  sets <- unlist(unname(collections), recursive = FALSE)
  sets <- sets[map_int(sets, \(s) sum(s %in% rownames(m))) >= min_size]
  multiScore(rankGenes(m), upSetColc = sets)$Scores
}
scores <- list(
  full = rbind(
    score(quant$E_norm$full, gs$full["GO:BP"], 15),
    score(quant$E_norm$full, c(gs$complexes, gs$clusters), 5)
  ),
  mito = score(quant$E_norm$mito, gs$mito, 10)
)
stopifnot(map_lgl(scores, \(s) identical(colnames(s), targets$sample_id)))

collection_of <- c(gs$full["GO:BP"], gs$mito, gs$complexes, gs$clusters) |>
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
    "Scores by mouse: GO:BP terms, CORUM complexes and STRING clusters.",
    "Scores by mouse: MitoCarta pathways in the mitochondrial universe.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

count(catalog, universe, collection)

sessionInfo()
