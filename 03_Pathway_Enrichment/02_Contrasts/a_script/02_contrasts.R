#!/usr/bin/env Rscript
# Which pathways each contrast moved: enrichVolcano runs fgsea (NES, the lead result), camera and
# fry on the limpa study, BH within each collection; full proteome against Hallmark, Reactome, KEGG
# and GO:BP, the mitochondrial universe against MitoCarta pathways. Redundant terms are flagged
# for display, never dropped.

source(here::here("R", "helpers.R"))
library(enrichVolcano)
BiocParallel::register(BiocParallel::SerialParam())

out <- stage_paths("03_Pathway_Enrichment", "02_Contrasts")
study_dir <- here("02_Differential_Expression", "02_Contrasts", "c_data")
inputs <- c(
  gene_sets = here("03_Pathway_Enrichment", "00_Gene_Sets", "c_data", "gene_sets.rds"),
  study_full = file.path(study_dir, "study_full", "da_results.csv"),
  study_mito = file.path(study_dir, "study_mito", "da_results.csv")
)

gs <- readRDS(inputs[["gene_sets"]])
collections <- list(full = gs$full, mito = gs$mito)
# MitoCarta pathways are smaller than MSigDB sets, so the mitochondrial run admits sets of 10.
min_size <- c(full = 15, mito = 10)

enrichment <- imap(collections, \(sets, universe) {
  study <- suppressMessages(read_study(file.path(study_dir, paste0("study_", universe)),
    species = "Mus musculus"
  ))
  set.seed(SEED)
  en <- suppressMessages(run_enrichment(study, sets,
    tests = c("fgsea", "camera", "fry"), inter_gene_cor = 0.01,
    min_size = min_size[[universe]], max_size = 500
  ))
  flat <- unlist(unname(sets), recursive = FALSE)
  map(en, \(e) dedup_terms(e, flat, method = "enrichmentmap"))
})

results <- imap(enrichment, \(tests, universe) {
  imap(tests, \(e, test) mutate(e@results, test = test, universe = universe, .before = 1))
}) |>
  unlist(recursive = FALSE) |>
  map(\(r) mutate(r, leading_edge = map_chr(leading_edge, \(x) paste(x, collapse = ";")))) |>
  list_rbind() |>
  left_join(select(gs$goslim, term, go_slim = slim), by = "term") |>
  mutate(ros_set = term %in% gs$ros) |>
  select(
    test, universe,
    collection = database, contrast, term, score, p, padj, size, direction,
    dedup_status, go_slim, ros_set, leading_edge
  )

summary <- results |>
  reframe(
    tier_counts(p, padj),
    share_true_effects = if (n() >= 100) 1 - limma::propTrueNull(p) else NA_real_,
    .by = c(test, universe, collection, contrast)
  )

nes_page <- function(x, title) {
  x <- slice_min(x, p, n = 15, with_ties = FALSE) |> arrange(score)
  par(mar = c(4, 22, 3, 1))
  barplot(x$score,
    names.arg = str_trunc(x$term, 55), horiz = TRUE, las = 1, cex.names = 0.55,
    col = if_else(x$padj <= ALPHA, "#B2182B", "grey70"), xlab = "NES", main = title,
    cex.main = 0.85
  )
  abline(v = 0)
}
plot_index <- map(CONTRASTS, \(cn) {
  fg <- results |>
    filter(test == "fgsea", contrast == cn, is.na(dedup_status) | dedup_status != "redundant")
  pages <- distinct(fg, universe, collection) |>
    pmap(\(universe, collection) {
      x <- filter(fg, universe == !!universe, collection == !!collection)
      \() nes_page(x, str_glue("{cn}: {collection}, 15 smallest p (red: FDR ≤ 0.05)"))
    })
  names(pages) <- str_glue("fgsea NES, {distinct(fg, collection)$collection}")
  write_pdf(file.path(out$reports, str_glue("02_contrasts_{cn}.pdf")), pages) |>
    mutate(contrast = cn, .before = 1)
}) |>
  list_rbind()

saveRDS(enrichment, file.path(out$data, "enrichment.rds"))
write_workbook(
  list(
    summary = summary, fgsea = filter(results, test == "fgsea"),
    camera = filter(results, test == "camera"), fry = filter(results, test == "fry"),
    ros_sets = filter(results, ros_set), plot_index = plot_index
  ),
  c(
    "Per test, universe, collection, contrast: p ≤ 0.05, chance, FDR tiers, true-effect share.",
    "fgsea: NES, p, BH FDR within collection, leading edge, redundancy flag, GO Slim ancestors.",
    "camera: competitive test with precision weights, inter-gene correlation 0.01.",
    "fry: self-contained rotation test with precision weights.",
    "The ROS sets named in advance, every test and contrast.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

summary |> filter(test == "fgsea", contrast %in% PRIMARY)

sessionInfo()
