# F03 a: pathway sets in the full proteome (GO:BP, Reactome, Hallmark, KEGG legacy) at fgsea
# FDR <= 0.05 in any contrast that collapsePathways keeps, grouped under the contrast that first
# selects them and drawn in all three (02_Enrichment).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"))
enr <- readRDS(inputs[["enrichment"]])

matrix <- contrast_pathway_matrix(
  filter(enr$results, universe == "full"),
  title = "Pathways in the full proteome",
  subtitle = "fgsea, FDR ≤ 0.05 in any contrast, redundant sets removed (collapsePathways)",
  nes_limit = 2.5
)

save_panel(matrix$plot, "F03_Pathways", "a_full_matrix", 110, 150)
list(
  plot = matrix$plot,
  data = left_join(matrix$cells, select(matrix$rows, term, collection, row_group), by = "term"),
  inputs = inputs
)
