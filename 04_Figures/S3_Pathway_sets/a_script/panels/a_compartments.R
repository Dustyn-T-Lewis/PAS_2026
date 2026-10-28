# S3 a: the GO Slim cellular compartments, fgsea per contrast in the full proteome, BH across the
# compartments (02_Enrichment), in F03's grid; rows ordered by the female NES, FDR <= 0.05 ringed.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"))
compartments <- readRDS(inputs[["enrichment"]])$compartments
rows <- compartments |>
  filter(contrast == "PAS_vs_VEH_female") |>
  arrange(desc(NES)) |>
  transmute(
    block = factor("GO Slim"), term = set, label = set_label(set, 70),
    shade = scales::alpha("grey94", 0.5)
  )

grid <- pathway_grid(
  transmute(compartments,
    term = set, column = contrast, value = NES, p = padj, selected = padj <= ALPHA
  ),
  rows,
  tibble(
    column = COMPARISON_CONTRASTS, column_block = "Contrast",
    label = CONTRAST_SHORT[COMPARISON_CONTRASTS]
  ),
  title = "Cellular compartments", subtitle = "GO Slim compartments, fgsea, full proteome",
  tones = c(Contrast = "Contrast", `GO Slim` = "GO Slim"),
  x_angle = 90, selected_label = paste("FDR ≤", ALPHA), fill_limit = 2.5, outline_shape = 21
)

height_mm <- 3 * nrow(rows) + 30
save_panel(grid, "S3_Pathway_sets", "a_compartments", 70, height_mm)
list(
  plot = grid,
  data = compartments |>
    left_join(select(rows, set = term, label), by = "set") |>
    transmute(set, label, contrast, NES, p, fdr = padj, size),
  inputs = inputs, height_mm = height_mm
)
