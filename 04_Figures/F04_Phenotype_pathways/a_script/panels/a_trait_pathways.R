# F04 a: the pathway sets that follow each phenotype within each sex: fgsea on each protein's
# association t with the trait under the within-treatment model (02_Enrichment), and the sets
# each trait shows (FDR <= 0.05, one per group of overlapping sets, both universes together). The
# two sexes share one set of rows in trait blocks: the ROWS_PER_SEX sets the trait shows in
# females, then those it shows in males, each by FDR then p, a set listed once in the figure. Each
# row is drawn against all six traits in both sexes, shaded across both grids for the sex with
# the smaller FDR at the block's trait, and its trait cell outlined in each sex that shows it.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"))
enr <- readRDS(inputs[["enrichment"]])
ROWS_PER_SEX <- 4
SEX_BLOCKS <- c(female = "Female traits", male = "Male traits")
# Each sex's block in its PAS contrast's tones; trait blocks fall back to grey boxes.
TONES <- c(`Female traits` = "PAS_vs_VEH_female", `Male traits` = "PAS_vs_VEH_male")
SEX_SHADES <- set_names(scales::alpha(SEX_COLOURS, 0.14), names(SEX_BLOCKS))

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
    label = label_with_database(set_label(term, 55), collection),
    shade = unname(SEX_SHADES[stronger_sex])
  )

drawn <- assoc |>
  inner_join(select(rows, term, row_block = block, stronger_sex), by = "term") |>
  mutate(selected = block == row_block & shown) |>
  arrange(sex, row_block, match(term, rows$term), block)
nes_limit <- ceiling(max(abs(drawn$NES)) * 10) / 10

grids <- imap(SEX_BLOCKS, \(column_block, s) {
  pathway_grid(
    drawn |>
      filter(sex == s) |>
      transmute(term, column = phenotype, NES, fdr = fgsea_fdr, selected),
    rows,
    tibble(column = names(TRAIT_SHORT), column_block, label = unname(TRAIT_SHORT)),
    tones = TONES, nes_limit = nes_limit, x_angle = 90, column_mm = 5,
    row_labels = s == "female", selected_label = "Chosen at this trait"
  )
})
plot <- wrap_elements(full = wrap_plots(grids, nrow = 1) + plot_layout(guides = "collect")) +
  labs(
    title = "Pathways per phenotype",
    subtitle = sprintf(
      paste(
        "Within-treatment model, fgsea FDR ≤ %s; one set per group of overlapping",
        "sets (most specific within %d-fold of the best p); up to %d per sex per trait, females",
        "first.\nRow shade: the sex with the smaller FDR at the block's trait"
      ),
      ALPHA, SPECIFIC_WINDOW, ROWS_PER_SEX
    )
  ) +
  theme_figure()

height_mm <- 3 * nrow(rows) + 36
save_panel(plot, "F04_Phenotype_pathways", "a_trait_pathways", 180, height_mm)
list(
  plot = plot,
  data = transmute(drawn,
    sex,
    block = as.character(row_block), term, collection, universe, phenotype, NES, fgsea_p,
    fgsea_fdr, shown, selected, stronger_sex, camera_p, camera_fdr
  ),
  inputs = inputs, height_mm = height_mm
)
