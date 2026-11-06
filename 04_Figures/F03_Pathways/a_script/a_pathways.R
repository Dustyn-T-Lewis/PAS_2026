# F03 a: the pathway sets each contrast shows (02_Enrichment's shown flag: fgsea FDR <= 0.05, one
# set per group of overlapping sets), each drawn in all three contrasts and at every trait in both
# sexes, the columns F04 draws. Full-proteome sets sit under the contrast that shows them, by FDR
# then p; MitoCarta sets (mitochondrial universe) shown in any contrast form the last block, by
# their smallest p. A set is listed once, in the first block that reaches it, and a block lists at
# most ROWS_PER_BLOCK sets. The contrast cells that put a row in its block are outlined: the
# block's contrast, or for MitoCarta every contrast showing the set; each row is shaded in its
# block's colour. Bars to the right count each contrast's distinct significant sets (FDR <= 0.05
# and shown) by direction and database: one pair per contrast block, and one pair per contrast in
# the MitoCarta block.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"))
enr <- readRDS(inputs[["enrichment"]])
ROWS_PER_BLOCK <- 6
BLOCKS <- c(CONTRAST_SHORT, MitoCarta = "MitoCarta")
# Each strip named for display, keyed to its contrast_strips() tones.
TONES <- c(set_names(names(BLOCKS), BLOCKS), SEX_TRAIT_TONES)
# Each row shaded in its block's strip fill, light enough that the dots stay readable.
SHADES <- scales::alpha(c(contrast_tone("fill"), MitoCarta = "grey85"), 0.5)
BAR_LABELS <- c(PAS_vs_VEH_female = "F", PAS_vs_VEH_male = "M", treatment_by_sex = "T×S")
WIDTH_MM <- 180

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
    label = label_with_database(set_label(term, 32), collection),
    shade = unname(SHADES[TONES[as.character(block)]])
  )

drawn <- contrast_trait_cells(enr, rows$term) |>
  inner_join(select(rows, term, block), by = "term") |>
  arrange(block, match(term, rows$term)) |>
  mutate(
    selected = shown & !is.na(contrast) &
      (block == CONTRAST_SHORT[contrast] | block == "MitoCarta")
  )
grid <- pathway_grid(
  select(drawn, term, column, value = NES, p = fdr, selected), rows,
  distinct(drawn, column, column_block, label),
  title = "Pathways per contrast",
  subtitle = sprintf(
    paste(
      "Rows chosen by contrast: fgsea FDR ≤ %s, one set per group of overlapping sets, up to %d",
      "per block. Traits: within-treatment association"
    ),
    ALPHA, ROWS_PER_BLOCK
  ),
  tones = TONES, selected_label = "Chosen in this contrast"
)

# Significant sets per contrast, direction and database: all of them, and the distinct ones the
# bars draw. Full-proteome sets count in their contrast's block, MitoCarta sets in the last.
counts <- enr$results |>
  filter(contrast %in% COMPARISON_CONTRASTS, padj <= ALPHA) |>
  mutate(block = if_else(universe == "mito", "MitoCarta", CONTRAST_SHORT[contrast])) |>
  summarise(
    distinct_fdr_le_0.05 = sum(shown), all_fdr_le_0.05 = n(),
    .by = c(block, contrast, direction, collection)
  ) |>
  arrange(factor(block, BLOCKS), match(contrast, COMPARISON_CONTRASTS), direction, collection)
sizes <- table(rows$block)
pairs <- tibble(
  contrast = rep(COMPARISON_CONTRASTS, 2),
  block = c(unname(CONTRAST_SHORT[COMPARISON_CONTRASTS]), rep("MitoCarta", 3))
) |>
  mutate(
    bar = BAR_LABELS[contrast], tone = contrast_tone("ink")[contrast],
    slot = if_else(block == "MitoCarta", match(contrast, COMPARISON_CONTRASTS), 1L),
    span = as.numeric(sizes[block]) / if_else(block == "MitoCarta", 3, 1),
    y = (slot - 0.5) * span + 0.5
  )
bars <- count_bars(
  transmute(counts, block,
    bar = BAR_LABELS[contrast], direction, collection, n = distinct_fdr_le_0.05
  ),
  pairs, rows
)

plot <- (grid + free(bars, type = "space", side = "b")) + plot_layout(guides = "collect")
height_mm <- grid_body_mm(rows) + 33
list(
  plot = plot,
  data = list(
    cells = transmute(drawn,
      block = as.character(block), set = term, collection, universe, contrast, sex, phenotype,
      NES, fdr, shown, selected
    ),
    counts = counts
  ),
  width_mm = WIDTH_MM, height_mm = height_mm
)
