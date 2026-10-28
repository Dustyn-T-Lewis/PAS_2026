# F03 a: the pathway sets each contrast shows (02_Enrichment's shown flag: fgsea FDR <= 0.05, one
# set per group of overlapping sets), each drawn in all three contrasts and at every trait in both
# sexes, the columns F04 draws. Full-proteome sets sit under the contrast that shows them, by FDR
# then p; MitoCarta sets (mitochondrial universe) shown in any contrast form the last block, by
# their smallest p. A set is listed once, in the first block that reaches it, and a block lists at
# most ROWS_PER_BLOCK sets. The contrast cells that put a row in its block are outlined: the
# block's contrast, or for MitoCarta every contrast showing the set; each row is shaded in its
# block's colour.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"))
enr <- readRDS(inputs[["enrichment"]])
ROWS_PER_BLOCK <- 6
BLOCKS <- c(CONTRAST_SHORT, MitoCarta = "MitoCarta")
# Each strip named for display, keyed to its contrast_strips() tones.
TONES <- c(set_names(names(BLOCKS), BLOCKS), SEX_TRAIT_TONES)
# Each row shaded in its block's strip fill, light enough that the dots stay readable.
SHADES <- scales::alpha(c(contrast_tone("fill"), MitoCarta = "grey85"), 0.5)

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
  mutate(block = factor(block, BLOCKS)) |>
  list_once(ROWS_PER_BLOCK) |>
  mutate(
    label = label_with_database(set_label(term, 60), collection),
    shade = unname(SHADES[TONES[as.character(block)]])
  )

drawn <- contrast_trait_cells(enr, rows$term) |>
  inner_join(select(rows, term, block), by = "term") |>
  arrange(block, match(term, rows$term)) |>
  mutate(
    selected = shown & !is.na(contrast) &
      (block == CONTRAST_SHORT[contrast] | block == "MitoCarta")
  )
plot <- pathway_grid(
  select(drawn, term, column, NES, fdr, selected), rows,
  distinct(drawn, column, column_block, label),
  title = "Pathways per contrast",
  subtitle = sprintf(
    paste(
      "Rows chosen by contrast: fgsea FDR ≤ %s, one set per group of overlapping sets, up to %d",
      "per block. Traits: within-treatment association"
    ),
    ALPHA, ROWS_PER_BLOCK
  ),
  tones = TONES, x_angle = 90, selected_label = "Chosen in this contrast"
)

height_mm <- 3 * nrow(rows) + nlevels(rows$block) + 30
save_panel(plot, "F03_Pathways", "a_pathways", 180, height_mm)
list(
  plot = plot,
  data = transmute(drawn,
    block = as.character(block), term, collection, universe, column_block, contrast, sex,
    phenotype, NES, p, fdr, shown, selected, camera_p, camera_fdr
  ),
  inputs = inputs, height_mm = height_mm
)
