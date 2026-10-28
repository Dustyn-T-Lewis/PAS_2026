# F04 a: the pathway sets that follow each phenotype within each sex: fgsea on each protein's
# association t with the trait under the within-treatment model (02_Enrichment), and the sets
# each trait shows (FDR <= 0.05, one per group of overlapping sets, both universes together). The
# rows sit in trait blocks: the ROWS_PER_SEX sets the trait shows in females, then those it shows
# in males, each by FDR then p, a set listed once in the figure. Each row is drawn in the three
# contrasts and at all six traits in both sexes, the columns F03 draws, shaded for the sex with
# the smaller FDR at the block's trait, and its trait cell outlined in each sex that shows it.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"))
enr <- readRDS(inputs[["enrichment"]])
ROWS_PER_SEX <- 4
SEX_SHADES <- set_names(scales::alpha(SEX_COLOURS, 0.14), c("female", "male"))

assoc <- filter(enr$association_fgsea, model == "treatment") |>
  mutate(block = factor(TRAIT_SHORT[phenotype], TRAIT_SHORT))
# Per set and trait, the sex with the smaller FDR (then p).
stronger <- assoc |>
  arrange(fgsea_fdr, fgsea_p) |>
  summarise(stronger_sex = first(sex), .by = c(block, term))
# Female then male candidates per trait, as list_once() blocks, so each sex takes its own quota.
rows <- assoc |>
  filter(shown) |>
  arrange(block, sex, fgsea_fdr, fgsea_p) |>
  mutate(trait = block, block = factor(paste(block, sex), unique(paste(block, sex)))) |>
  list_once(ROWS_PER_SEX) |>
  transmute(block = trait, term, collection) |>
  left_join(stronger, by = c("block", "term")) |>
  mutate(
    label = label_with_database(set_label(term, 60), collection),
    shade = unname(SEX_SHADES[stronger_sex])
  )

drawn <- contrast_trait_cells(enr, rows$term) |>
  inner_join(select(rows, term, block, stronger_sex), by = "term") |>
  arrange(block, match(term, rows$term)) |>
  mutate(selected = shown & !is.na(phenotype) & block == TRAIT_SHORT[phenotype])
plot <- pathway_grid(
  select(drawn, term, column, NES, fdr, selected), rows,
  distinct(drawn, column, column_block, label),
  title = "Pathways per phenotype",
  subtitle = sprintf(
    paste(
      "Rows chosen by trait: within-treatment fgsea FDR ≤ %s, one set per group of overlapping",
      "sets, up to %d per sex per trait. Row shade: the sex with the smaller FDR there"
    ),
    ALPHA, ROWS_PER_SEX
  ),
  tones = SEX_TRAIT_TONES, x_angle = 90, selected_label = "Chosen at this trait"
)

height_mm <- 3 * nrow(rows) + nlevels(rows$block) + 30
save_panel(plot, "F04_Phenotype_pathways", "a_trait_pathways", 180, height_mm)
list(
  plot = plot,
  data = transmute(drawn,
    block = as.character(block), term, collection, universe, column_block, contrast, sex,
    phenotype, NES, p, fdr, shown, selected, stronger_sex, camera_p, camera_fdr
  ),
  inputs = inputs, height_mm = height_mm
)
