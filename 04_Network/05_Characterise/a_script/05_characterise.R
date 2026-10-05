#!/usr/bin/env Rscript
# What each module is: over-represented gene sets (fgsea::fora against the detected proteome, BH
# within collection), its hub proteins, and whether its members interact more than expected in
# STRING (STRINGdb PPI enrichment against the detected proteome). STRING's answers are cached.

source(here::here("R", "helpers.R"))

out <- stage_paths("04_Network", "05_Characterise")
study_dir <- here("02_Differential_Expression", "02_Contrasts", "c_data", "study_full")
inputs <- c(
  modules = here("04_Network", "01_Build_Modules", "c_data", "modules.rds"),
  gene_sets = here("03_Pathway_Enrichment", "00_Gene_Sets", "c_data", "gene_sets.rds"),
  study = file.path(study_dir, "da_results.csv")
)
STRING_FILES <- here("04_Network", "stringdb_cache")
CACHE <- file.path(out$data, "string_cache.rds")

mods <- readRDS(inputs[["modules"]])
gs <- readRDS(inputs[["gene_sets"]])
mapping <- suppressMessages(enrichVolcano::read_study(study_dir, species = "Mus musculus"))$da |>
  distinct(protein, gene)
members <- mods$modules |>
  filter(module != "grey") |>
  left_join(mapping, by = "protein", suffix = c("_diann", "")) |>
  split(~module)
universe <- unique(mapping$gene)

annotation <- imap(members, \(m, module) {
  imap(gs$full, \(sets, collection) {
    fgsea::fora(sets, unique(m$gene), universe, minSize = 15, maxSize = 500) |>
      as_tibble() |>
      transmute(
        module = module, collection = collection, term = pathway, overlap, size, p = pval,
        padj, overlap_genes = map_chr(overlapGenes, \(g) paste(g, collapse = ";"))
      )
  }) |>
    list_rbind()
}) |>
  list_rbind() |>
  arrange(module, p)

key <- tools::md5sum(inputs[["modules"]])
if (!file.exists(CACHE) || !identical(readRDS(CACHE)$key, unname(key))) {
  dir.create(STRING_FILES, showWarnings = FALSE)
  string_db <- STRINGdb::STRINGdb$new(
    version = "12.0", species = 10090, score_threshold = 400,
    network_type = "full", input_directory = STRING_FILES
  )
  string_map <- string_db$map(as.data.frame(mods$modules["protein"]), "protein",
    removeUnmappedRows = TRUE
  )
  string_db$set_background(unique(string_map$STRING_id))
  string <- imap(members, \(m, module) {
    ids <- unique(string_map$STRING_id[string_map$protein %in% m$protein])
    e <- string_db$get_ppi_enrichment(ids)
    tibble(
      module = module, mapped = length(ids), edges = e$edges, expected = e$lambda,
      p = e$enrichment
    )
  }) |>
    list_rbind()
  saveRDS(list(key = unname(key), string = string, fetched = Sys.Date()), CACHE)
}
string <- readRDS(CACHE)$string |>
  mutate(fold = edges / expected, padj = p.adjust(p, "BH"))

top_terms <- annotation |>
  slice_min(p, n = 3, by = module, with_ties = FALSE) |>
  summarise(top_terms = paste(term, collapse = "; "), .by = module)
overview <- mods$sizes |>
  filter(module != "grey") |>
  left_join(top_terms, by = "module") |>
  left_join(select(string, module, string_fold = fold, string_padj = padj), by = "module") |>
  left_join(summarise(mods$hubs, hubs = paste(head(gene, 5), collapse = ", "), .by = module),
    by = "module"
  )

plot_index <- write_pdf(file.path(out$reports, "05_characterise.pdf"), list(
  "STRING interaction enrichment per module" = \() {
    par(mar = c(4, 8, 3, 1))
    x <- arrange(string, fold)
    barplot(log2(x$fold),
      names.arg = x$module, horiz = TRUE, las = 1,
      col = if_else(x$padj <= ALPHA, "#B2182B", "grey70"), xlab = "log2 observed / expected edges",
      main = "STRING PPI enrichment, detected-proteome background (red: FDR ≤ 0.05)"
    )
    abline(v = 0)
  }
))

write_workbook(
  list(
    overview = overview, annotation = annotation, hubs = mods$hubs, string = string,
    plot_index = plot_index
  ),
  c(
    "Per module: size, MitoCarta share, top three terms, STRING fold and FDR, top five hubs.",
    "Over-represented sets per module (fgsea::fora, detected background, BH within collection).",
    "The ten highest-kME proteins of each module.",
    "STRING PPI enrichment per module: mapped proteins, edges, expected, fold, p, BH FDR.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

overview

sessionInfo()
