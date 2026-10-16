# F03 a: the pathway sets each contrast shows (02_Enrichment's shown flag: fgsea FDR <= 0.05, one
# set per group of overlapping sets), each drawn in all three contrasts. Full-proteome sets sit
# under the contrast that shows them, by FDR then p; MitoCarta sets (mitochondrial universe) shown
# in any contrast form the last block, by their smallest p. A set is listed once, in the first
# block that reaches it, and a block lists at most ROWS_PER_BLOCK sets.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"))
enr <- readRDS(inputs[["enrichment"]])
ROWS_PER_BLOCK <- 10
# Each block and column strip named for display, keyed to its contrast_strips() tones.
TONES <- c(set_names(names(CONTRAST_SHORT), CONTRAST_SHORT), MitoCarta = "MitoCarta")

by_contrast <- enr$results |>
  filter(universe == "full", shown) |>
  arrange(factor(contrast, COMPARISON_CONTRASTS), padj, p) |>
  transmute(block = CONTRAST_SHORT[contrast], term, collection)
mitocarta <- enr$results |>
  filter(universe == "mito", any(shown), .by = term) |>
  summarise(p = min(p), .by = c(term, collection)) |>
  arrange(p) |>
  transmute(block = "MitoCarta", term, collection)
rows <- bind_rows(by_contrast, mitocarta) |>
  mutate(block = factor(block, names(TONES))) |>
  list_once(ROWS_PER_BLOCK) |>
  mutate(label = set_label(term, 60))

drawn <- enr$results |>
  inner_join(select(rows, term, block), by = "term") |>
  arrange(block, match(term, rows$term), factor(contrast, COMPARISON_CONTRASTS))
plot <- pathway_grid(
  transmute(drawn, term, column = contrast, NES, fdr = padj), rows,
  tibble(column = COMPARISON_CONTRASTS, column_block = CONTRAST_SHORT, label = ""),
  title = "Pathways per contrast",
  subtitle = sprintf(
    paste(
      "fgsea FDR ≤ %s; one set per group of overlapping sets (most specific within %d-fold of",
      "the best p); up to %d per block"
    ),
    ALPHA, SPECIFIC_WINDOW, ROWS_PER_BLOCK
  ),
  tones = TONES, column_mm = 26
) +
  theme(axis.text.x = element_blank())

height_mm <- 3 * nrow(rows) + nlevels(rows$block) + 18
save_panel(plot, "F03_Pathways", "a_pathways", 180, height_mm)
list(
  plot = plot,
  data = transmute(drawn,
    block = as.character(block), term, collection, universe, contrast, NES, p, padj, shown,
    camera_p, camera_fdr
  ),
  inputs = inputs, height_mm = height_mm
)
