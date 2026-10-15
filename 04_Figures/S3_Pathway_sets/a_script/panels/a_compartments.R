# S3 a: the GO Slim cellular compartments, fgsea per contrast in the full proteome, BH across the
# compartments (02_Enrichment); rows ordered by the female NES.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"))
compartments <- readRDS(inputs[["enrichment"]])$compartments
order <- compartments |>
  filter(contrast == "PAS_vs_VEH_female") |>
  arrange(NES) |>
  pull(set)
cells <- compartments |>
  transmute(
    contrast, NES,
    fdr = padj, size,
    row = factor(set, order),
    label = str_to_sentence(str_replace_all(str_remove(set, "^GOSLIM_"), "_", " ")),
    row_group = ""
  )

plot <- set_dot_matrix(cells,
  title = "Cellular compartments",
  subtitle = "GO Slim compartments, fgsea, full proteome"
)

save_panel(plot, "S3_Pathway_sets", "a_compartments", 80, 120)
list(plot = plot, data = select(cells, -row), inputs = inputs)
