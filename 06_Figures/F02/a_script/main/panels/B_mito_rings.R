#!/usr/bin/env Rscript
# F02 b. enrichVolcano rings for the mitochondrial proteome, female above male: the 12 strongest
# collapsed fgsea terms each (03_Pathway_Enrichment/08_Display), proteins at FDR ≤ 0.10 coloured.

source(here::here("R", "panels.R"))

SEXES <- c(PAS_vs_VEH_female = "Female", PAS_vs_VEH_male = "Male")

inputs <- c(
  display = here("03_Pathway_Enrichment", "08_Display", "c_data", "display.rds"),
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
display <- readRDS(inputs[["display"]])
hits <- filter(readRDS(inputs[["contrasts"]])$results, universe == "mito", adj.P.Val <= 0.10)
da <- study_da("mito")
rings <- imap(SEXES, \(sex, cn) {
  terms <- filter(display$ring_terms, universe == "mito", contrast == cn)
  sex_hits <- filter(hits, contrast == cn)
  pathway_ring(display, da, "mito", cn,
    title = str_glue("{sex}, mitochondrial proteome"),
    subtitle = str_glue(
      "{RING_POOLS[['mito']]}\n",
      "{nrow(terms)} non-redundant terms, {sum(terms$padj <= ALPHA)} at FDR ≤ 0.05 (*)\n",
      "Proteins at FDR ≤ 0.10: {sum(sex_hits$logFC < 0)} down, {sum(sex_hits$logFC > 0)} up"
    )
  )
})
panel <- rings[[1]] / rings[[2]]
panel_data <- filter(display$ring_terms, universe == "mito", contrast %in% names(SEXES)) |>
  select(contrast, database, term, NES = score, p, padj)
panel_note <- "Ring terms per sex, mitochondrial proteome (03_Pathway_Enrichment/08_Display)."
save_panel(panel, "F02", "main", "B_mito_rings", 100, 190)
