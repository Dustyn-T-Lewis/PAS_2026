# F03 b: MitoCarta pathways in the mitochondrial universe, selected and drawn as panel a.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"))
enr <- readRDS(inputs[["enrichment"]])

matrix <- contrast_pathway_matrix(
  filter(enr$results, universe == "mito"),
  title = "MitoCarta pathways",
  subtitle = "Mitochondrial universe; selected as in a",
  nes_limit = 2.5
)

save_panel(matrix$plot, "F03_Pathways", "b_mitocarta_matrix", 70, 55)
list(
  plot = matrix$plot,
  data = left_join(matrix$cells, select(matrix$rows, term, collection, row_group), by = "term"),
  inputs = inputs
)
