# F04 a: the pathway sets that follow each phenotype within each sex: fgsea on each protein's
# association t with the trait under the within-treatment model (02_Enrichment), and the sets
# each trait shows (FDR <= 0.05, one per group of overlapping sets, both universes together). The
# rows sit in trait blocks: the ROWS_PER_SEX sets the trait shows in females, then those it shows
# in males, each by FDR then p, a set listed once in the figure. Each row is drawn in the three
# contrasts and at all six traits in both sexes, the columns F03 draws, shaded for the sex that
# chose it, and its trait cell outlined in each sex that shows it. Bars to the right count each
# trait's distinct significant sets (FDR <= 0.05 and shown) by direction and database, the female
# pair beside the female rows and the male pair beside the male rows.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"))
enr <- readRDS(inputs[["enrichment"]])
ROWS_PER_SEX <- 4
SEX_SHADES <- set_names(scales::alpha(SEX_COLOURS, 0.14), c("female", "male"))
WIDTH_MM <- 212

assoc <- filter(enr$association_fgsea, model == "treatment") |>
  mutate(block = factor(TRAIT_SHORT[phenotype], TRAIT_SHORT))
# Female then male candidates per trait, as list_once() blocks, so each sex takes its own quota.
rows <- assoc |>
  filter(shown) |>
  arrange(block, sex, fgsea_fdr, fgsea_p) |>
  mutate(trait = block, block = factor(paste(block, sex), unique(paste(block, sex)))) |>
  list_once(ROWS_PER_SEX) |>
  transmute(
    block = trait, term, collection, sex,
    label = label_with_database(set_label(term, 60), collection),
    shade = unname(SEX_SHADES[sex])
  )

drawn <- contrast_trait_cells(enr, rows$term) |>
  inner_join(select(rows, term, block, chosen_sex = sex), by = "term") |>
  arrange(block, match(term, rows$term)) |>
  mutate(selected = shown & !is.na(phenotype) & block == TRAIT_SHORT[phenotype])
grid <- pathway_grid(
  select(drawn, term, column, NES, fdr, selected), rows,
  distinct(drawn, column, column_block, label),
  title = "Pathways per phenotype",
  subtitle = sprintf(
    paste(
      "Rows chosen by trait: within-treatment fgsea FDR ≤ %s, one set per group of overlapping",
      "sets, up to %d per sex per trait, shaded for the sex that chose them"
    ),
    ALPHA, ROWS_PER_SEX
  ),
  tones = SEX_TRAIT_TONES, x_angle = 90, selected_label = "Chosen at this trait"
)

# Significant sets per trait, sex, direction and database: all of them, and the distinct ones
# the bars draw.
counts <- assoc |>
  filter(fgsea_fdr <= ALPHA) |>
  mutate(direction = if_else(NES > 0, "up", "down")) |>
  summarise(
    n_distinct = sum(shown), n_all = n(), .by = c(block, phenotype, sex, direction, collection)
  ) |>
  arrange(block, sex, direction, collection) |>
  mutate(block = as.character(block))
# Each sex's pair centred on the rows that sex chose; a sex with none takes the block's last row.
block_rows <- table(rows$block)
pairs <- rows |>
  mutate(row = row_number(), .by = block) |>
  summarise(y = mean(row), span = n(), .by = c(block, sex)) |>
  mutate(block = as.character(block)) |>
  complete(block, sex = c("female", "male"), fill = list(span = 1L)) |>
  mutate(
    y = coalesce(y, as.numeric(block_rows[block])), bar = if_else(sex == "female", "F", "M"),
    tone = unname(SEX_COLOURS[str_to_title(sex)])
  )
bars <- count_bars(
  transmute(counts, block,
    bar = if_else(sex == "female", "F", "M"), direction, collection,
    n = n_distinct
  ),
  pairs, rows
)

plot <- (grid + free(bars, type = "space", side = "b")) + plot_layout(guides = "collect")
height_mm <- 3 * nrow(rows) + nlevels(rows$block) + 32
save_panel(plot, "F04_Phenotype_pathways", "a_trait_pathways", WIDTH_MM, height_mm)
list(
  plot = plot,
  data = transmute(drawn,
    block = as.character(block), term, collection, universe, column_block, contrast, sex,
    phenotype, NES, p, fdr, shown, selected, chosen_sex, camera_p, camera_fdr
  ),
  counts = select(counts, -block), inputs = inputs, width_mm = WIDTH_MM,
  height_mm = height_mm
)
