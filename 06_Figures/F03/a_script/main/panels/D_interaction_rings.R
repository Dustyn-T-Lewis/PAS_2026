#!/usr/bin/env Rscript
# F03 d. enrichVolcano rings for the treatment-by-sex contrast, full and mitochondrial proteome
# (03_Pathway_Enrichment/08_Display); > 0 means PAS raises the term more in females.

source(here::here("R", "panels.R"))

UNIVERSE_LABELS <- c(full = "full proteome", mito = "mitochondrial proteome")

inputs <- c(display = here("03_Pathway_Enrichment", "08_Display", "c_data", "display.rds"))
display <- readRDS(inputs[["display"]])
rings <- map(names(UNIVERSE_LABELS), \(u) {
  terms <- filter(display$ring_terms, universe == u, contrast == "treatment_by_sex")
  pathway_ring(display, study_da(u), u, "treatment_by_sex",
    title = str_glue("Treatment × sex, {UNIVERSE_LABELS[[u]]}"),
    subtitle = str_glue(
      "{nrow(terms)} non-redundant terms, {sum(terms$padj <= ALPHA)} at FDR ≤ 0.05 (*);\n",
      "> 0: PAS raises the term more in females than in males"
    )
  )
})
panel <- rings[[1]] / rings[[2]] + plot_layout(guides = "collect")
panel_data <- filter(display$ring_terms, contrast == "treatment_by_sex") |>
  select(universe, database, term, NES = score, p, padj)
panel_note <- "Interaction ring terms per universe (03_Pathway_Enrichment/08_Display)."
save_panel(panel, "F03", "main", "D_interaction_rings", 110, 190)
