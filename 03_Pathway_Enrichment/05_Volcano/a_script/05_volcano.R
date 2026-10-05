#!/usr/bin/env Rscript
# enrichVolcano rings: each contrast's protein volcano with its strongest fgsea terms around it,
# redundant terms left out by the EnrichmentMap flag. One PDF per contrast, a page per collection.

source(here::here("R", "helpers.R"))
library(enrichVolcano)

out <- stage_paths("03_Pathway_Enrichment", "05_Volcano")
study_dir <- here("02_Differential_Expression", "02_Contrasts", "c_data")
inputs <- c(
  enrichment = here("03_Pathway_Enrichment", "02_Contrasts", "c_data", "enrichment.rds"),
  study_full = file.path(study_dir, "study_full", "da_results.csv"),
  study_mito = file.path(study_dir, "study_mito", "da_results.csv")
)

enrichment <- readRDS(inputs[["enrichment"]])
da <- map(c(full = "full", mito = "mito"), \(u) {
  path <- file.path(study_dir, paste0("study_", u))
  suppressMessages(read_study(path, species = "Mus musculus"))$da
})
panels <- tribble(
  ~universe, ~collection,
  "full", "Hallmark", "full", "Reactome", "full", "KEGG", "full", "GO:BP", "mito", "MitoCarta"
)

plot_index <- map(CONTRASTS, \(cn) {
  pages <- pmap(panels, \(universe, collection) {
    p <- plot_volcano_ring(da[[universe]], enrichment[[universe]]$fgsea,
      contrast = cn, databases = collection, collapse = TRUE, term_threshold = ALPHA,
      n_terms = 12, title = str_glue("{cn}: {collection}"),
      subtitle = "fgsea terms at FDR ≤ 0.05, one per EnrichmentMap cluster"
    )
    \() print(p)
  })
  names(pages) <- str_glue("volcano ring, {panels$collection}")
  file <- file.path(out$reports, str_glue("05_volcano_{cn}.pdf"))
  write_pdf(file, pages, width = 8, height = 8) |>
    mutate(contrast = cn, .before = 1)
}) |>
  list_rbind()

write_workbook(
  list(plot_index = plot_index),
  "Which PDF page shows which contrast and collection.",
  out$workbook, inputs
)

sessionInfo()
