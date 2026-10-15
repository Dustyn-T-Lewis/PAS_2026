# F04 a: pathway sets tracking each phenotype across the ten female mice. fgsea on each protein's
# association t with the trait, pooled over treatment and within treatment (02_Enrichment); per
# trait, the four kept sets at FDR <= 0.05 with the smallest FDR under either model that no
# earlier trait already lists.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"))
enr <- readRDS(inputs[["enrichment"]])

matrix <- phenotype_pathway_matrix(enr$association_fgsea, "female",
  title = "Females: pathways tracking each phenotype",
  subtitle = "fgsea on protein–trait association; each set under the first trait to select it",
  nes_limit = 3
)

save_panel(matrix$plot, "F04_Phenotype_pathways", "a_female", 90, 150)
list(
  plot = matrix$plot,
  data = left_join(matrix$cells, select(matrix$rows, term, collection, row_group),
    by = "term", relationship = "many-to-many"
  ),
  inputs = inputs
)
