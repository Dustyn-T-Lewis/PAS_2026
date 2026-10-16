# F04 a: the pathway sets that follow each phenotype within each sex: fgsea on each protein's
# association t with the trait under the within-treatment model (02_Enrichment), and the sets
# each trait shows (FDR <= 0.05, one per group of overlapping sets, both universes together), by
# FDR then p. Each sex is its own grid of trait blocks, a set listed once per sex and at most
# ROWS_PER_TRAIT sets per trait, drawn against all six traits; the two grids share one key.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"))
enr <- readRDS(inputs[["enrichment"]])
ROWS_PER_TRAIT <- 5
SEX_BLOCKS <- c(female = "Female traits", male = "Male traits")
# Each sex's block in its PAS contrast's tones; trait blocks fall back to grey boxes.
TONES <- c(`Female traits` = "PAS_vs_VEH_female", `Male traits` = "PAS_vs_VEH_male")

assoc <- filter(enr$association_fgsea, model == "treatment")
rows <- map(set_names(names(SEX_BLOCKS)), \(s) {
  assoc |>
    filter(sex == s, shown) |>
    arrange(fgsea_fdr, fgsea_p) |>
    transmute(block = factor(TRAIT_SHORT[phenotype], TRAIT_SHORT), term, collection) |>
    list_once(ROWS_PER_TRAIT) |>
    mutate(label = set_label(term, 40))
})
drawn <- imap(rows, \(r, s) {
  assoc |>
    filter(sex == s) |>
    inner_join(select(r, term, block), by = "term") |>
    arrange(block, match(term, r$term), factor(phenotype, names(TRAIT_SHORT)))
})
nes_limit <- ceiling(max(abs(bind_rows(drawn)$NES)) * 10) / 10
# Row slots per trait block: the larger of the two sexes, so the blocks line up across.
slots <- reduce(map(rows, \(r) table(r$block)), pmax)

grids <- imap(rows, \(r, s) {
  pathway_grid(
    transmute(drawn[[s]], term, column = phenotype, NES, fdr = fgsea_fdr), r,
    tibble(
      column = names(TRAIT_SHORT), column_block = SEX_BLOCKS[[s]], label = unname(TRAIT_SHORT)
    ),
    tones = TONES, nes_limit = nes_limit, x_angle = 90, column_mm = 4.5, block_rows = slots
  )
})
plot <- wrap_elements(full = wrap_plots(grids, nrow = 1) + plot_layout(guides = "collect")) +
  labs(
    title = "Pathways that follow each phenotype",
    subtitle = sprintf(
      paste(
        "Within-treatment model, fgsea FDR ≤ %s; one set per group of overlapping sets (most",
        "specific within %d-fold of the best p); up to %d per trait"
      ),
      ALPHA, SPECIFIC_WINDOW, ROWS_PER_TRAIT
    )
  ) +
  theme_figure()

height_mm <- 3 * sum(slots) + 36
save_panel(plot, "F04_Phenotype_pathways", "a_trait_pathways", 180, height_mm)
list(
  plot = plot,
  data = drawn |>
    bind_rows() |>
    transmute(
      sex,
      block = as.character(block), term, collection, universe, phenotype, NES, fgsea_p,
      fgsea_fdr, shown, camera_p, camera_fdr
    ),
  inputs = inputs, height_mm = height_mm
)
