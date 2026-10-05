#!/usr/bin/env Rscript
# F02 a. enrichVolcano rings for the full proteome, female above male: the 12 strongest collapsed
# fgsea terms each (03_Pathway_Enrichment/08_Display), proteins at FDR ≤ 0.10 coloured.

source(here::here("R", "panels.R"))

SEXES <- c(PAS_vs_VEH_female = "Female", PAS_vs_VEH_male = "Male")

inputs <- c(
  display = here("03_Pathway_Enrichment", "08_Display", "c_data", "display.rds"),
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
display <- readRDS(inputs[["display"]])
hits <- filter(readRDS(inputs[["contrasts"]])$results, universe == "full", adj.P.Val <= 0.10)
da <- study_da("full")
rings <- imap(SEXES, \(sex, cn) {
  terms <- filter(display$ring_terms, universe == "full", contrast == cn)
  sex_hits <- filter(hits, contrast == cn)
  pathway_ring(display, da, "full", cn,
    title = str_glue("{sex}, full proteome"),
    subtitle = str_glue(
      "{RING_POOLS[['full']]}\n",
      "{nrow(terms)} non-redundant terms, {sum(terms$padj <= ALPHA)} at FDR ≤ 0.05 (*)\n",
      "Proteins at FDR ≤ 0.10: {sum(sex_hits$logFC < 0)} down, {sum(sex_hits$logFC > 0)} up"
    )
  )
})
panel <- rings[[1]] / rings[[2]]
panel_data <- filter(display$ring_terms, universe == "full", contrast %in% names(SEXES)) |>
  select(contrast, database, term, NES = score, p, padj)
panel_note <- "Ring terms per sex, full proteome (03_Pathway_Enrichment/08_Display)."
save_panel(panel, "F02", "main", "A_full_rings", 100, 190)
