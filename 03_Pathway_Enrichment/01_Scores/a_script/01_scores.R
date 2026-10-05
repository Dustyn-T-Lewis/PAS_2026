#!/usr/bin/env Rscript
# One singscore per mouse per gene set, from the same collections and the same one-protein-per-gene
# rule as the pathway tests: full-proteome collections on the full matrix, MitoCarta on the
# mitochondrial matrix. A score is the mean within-mouse rank of the set's members and tests
# nothing; association and classification steps read it.

source(here::here("R", "helpers.R"))
library(singscore)

out <- stage_paths("03_Pathway_Enrichment", "01_Scores")
inputs <- c(
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds"),
  gene_sets = here("03_Pathway_Enrichment", "00_Gene_Sets", "c_data", "gene_sets.rds")
)

de <- readRDS(inputs[["contrasts"]])
gs <- readRDS(inputs[["gene_sets"]])
# The symbols and abundances enrichVolcano itself uses, so a set scored here has the members it
# was tested with.
study_dir <- here("02_Differential_Expression", "02_Contrasts", "c_data", "study_full")
mapping <- suppressMessages(enrichVolcano::read_study(study_dir, species = "Mus musculus"))$da |>
  distinct(protein, gene, abundance)
symbols <- set_names(mapping$gene, mapping$protein)
abundance <- set_names(mapping$abundance, mapping$protein)

by_gene <- function(E) {
  keep <- one_per_gene(rownames(E), symbols[rownames(E)], abundance[rownames(E)])
  m <- E[keep, ]
  rownames(m) <- symbols[rownames(m)]
  m
}
score <- function(E, collections, min_size) {
  m <- by_gene(E)
  sets <- unlist(unname(collections), recursive = FALSE)
  sets <- sets[map_int(sets, \(s) sum(s %in% rownames(m))) >= min_size]
  multiScore(rankGenes(m), upSetColc = sets)$Scores
}
scores <- list(full = score(de$E_norm$full, gs$full, 15), mito = score(de$E_norm$mito, gs$mito, 10))
stopifnot(all(map_lgl(scores, \(s) identical(colnames(s), rownames(de$targets)))))

collection_of <- c(gs$full, gs$mito) |>
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
    m <- scores[[universe]][rows, ]
    title <- str_glue("{collection} singscores, z-scored by set, rows clustered (display only)")
    \() ComplexHeatmap::draw(overview_heatmap(m, de$targets$group, title))
  })
names(pages) <- str_glue("singscore overview, {distinct(catalog, universe, collection)$collection}")
plot_index <- write_pdf(file.path(out$reports, "01_scores_overview.pdf"), pages)

saveRDS(list(scores = scores, collection = collection_of), file.path(out$data, "set_scores.rds"))
write_workbook(
  list(
    catalog = catalog,
    scores_full = as_tibble(scores$full, rownames = "set"),
    scores_mito = as_tibble(scores$mito, rownames = "set"),
    plot_index = plot_index
  ),
  c(
    "Every scored set: universe, collection, mean score.",
    "Set by mouse scores, full-proteome collections.",
    "Set by mouse scores, MitoCarta pathways in the mitochondrial universe.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

count(catalog, universe, collection)

sessionInfo()
