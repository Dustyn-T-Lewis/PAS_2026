#!/usr/bin/env Rscript
# Which cell compartments each contrast moved: fgsea (competitive) and fry (self-contained: does the
# compartment's protein set move as a whole) on the GO Slim cellular-component sets over the full
# proteome, BH within the collection. The test is competitive, so a compartment rises or falls
# against the rest of the proteome: a shift in share, not an absolute amount. The pellet is
# mitochondria-enriched and its MitoCarta share differs by sex (see 02_Differential_Expression).

source(here::here("R", "helpers.R"))
BiocParallel::register(BiocParallel::SerialParam())

out <- stage_paths("03_Pathway_Enrichment", "07_Compartments")
inputs <- c(
  gene_sets = here("03_Pathway_Enrichment", "00_Gene_Sets", "c_data", "gene_sets.rds"),
  study_full = here(
    "02_Differential_Expression", "02_Contrasts", "c_data", "study_full", "da_results.csv"
  )
)
gs <- readRDS(inputs[["gene_sets"]])
study <- suppressMessages(enrichVolcano::read_study(dirname(inputs[["study_full"]]),
  species = "Mus musculus"
))

# No upper size limit: the large compartments (cytosol, nucleus) are the point of the test.
set.seed(SEED)
tests <- suppressMessages(enrichVolcano::run_enrichment(study, gs$compartments,
  tests = c("fgsea", "fry"), min_size = 15, max_size = Inf
))
compartments <- tests$fgsea
saveRDS(compartments, file.path(out$data, "compartments.rds"))
saveRDS(tests$fry, file.path(out$data, "compartments_fry.rds"))

results <- compartments@results |>
  select(contrast, term, size, NES = score, p, padj, leading_edge) |>
  mutate(leading_edge = map_chr(leading_edge, \(x) paste(x, collapse = ";")))
fry <- tests$fry@results |>
  select(contrast, term, size, direction, score, p, padj)
write_workbook(
  list(fgsea = results, fry = fry),
  c(
    "fgsea per contrast and GO Slim compartment: NES, p, BH FDR in the collection, leading edge.",
    "fry per contrast and compartment: direction, score, p, BH FDR within the collection."
  ),
  out$workbook, inputs
)

results |>
  filter(contrast %in% PRIMARY, padj <= ALPHA) |>
  arrange(contrast, padj)
sessionInfo()
