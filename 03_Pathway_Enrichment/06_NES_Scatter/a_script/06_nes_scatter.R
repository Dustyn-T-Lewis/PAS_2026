#!/usr/bin/env Rscript
# Pathway effects in one sex against the other: female PAS against male PAS, and the sex
# difference without against with treatment. enrichVolcano's NES scatter, a page per collection;
# whether the sexes differ is tested by treatment_by_sex in 02_Contrasts, not read off this plot.

source(here::here("R", "helpers.R"))
library(enrichVolcano)

out <- stage_paths("03_Pathway_Enrichment", "06_NES_Scatter")
inputs <- c(enrichment = here("03_Pathway_Enrichment", "02_Contrasts", "c_data", "enrichment.rds"))
enrichment <- readRDS(inputs[["enrichment"]])

PAIRS <- list(
  treatment = c(x = "PAS_vs_VEH_male", y = "PAS_vs_VEH_female"),
  sex = c(x = "sex_in_PAS", y = "sex_in_VEH")
)
panels <- tribble(
  ~universe, ~collection,
  "full", "Hallmark", "full", "Reactome", "full", "GO:BP", "mito", "MitoCarta"
)

plot_index <- imap(PAIRS, \(pair, name) {
  pages <- pmap(panels, \(universe, collection) {
    p <- plot_scatter(enrichment[[universe]]$fgsea,
      x = pair[["x"]], y = pair[["y"]], comparison = "concordance", databases = collection,
      collapse = TRUE, term_threshold = ALPHA, label_n = 10
    )
    \() print(p)
  })
  names(pages) <- str_glue("NES scatter, {panels$collection}")
  file <- file.path(out$reports, str_glue("06_nes_scatter_{name}.pdf"))
  write_pdf(file, pages, width = 8, height = 8) |>
    mutate(pair = str_glue("{pair[['y']]} against {pair[['x']]}"), .before = 1)
}) |>
  list_rbind()

write_workbook(
  list(plot_index = plot_index),
  "Which PDF page shows which pair of contrasts and collection.",
  out$workbook, inputs
)

sessionInfo()
