# S5 b: the literature-informed and ROS sets, grouped by the hypothesis that named them, with the
# direction each hypothesis expected; fgsea per contrast, BH across the list (02_Enrichment), in
# F03's grid with each set's database as a square, read from its MSigDB or MitoCarta prefix.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"))
literature <- readRDS(inputs[["enrichment"]])$literature
PREFIX_DATABASE <- c(
  GOBP = "GO:BP", HALLMARK = "Hallmark", MITOCARTA = "MitoCarta", REACTOME = "Reactome"
)
# Neighbouring blocks banded in alternating greys; the one-set block's row taller, so its
# two-line strip fits.
SHADES <- scales::alpha(c("grey94", "grey86"), 0.5)
MIN_BLOCK_MM <- 5.5

rows <- literature |>
  filter(contrast == "PAS_vs_VEH_female") |>
  mutate(block = factor(str_wrap(hypothesis, 12))) |>
  arrange(block) |>
  transmute(
    block,
    term = set, label = sprintf("%s [%s]", set_label(set, 70), expected),
    collection = unname(PREFIX_DATABASE[str_extract(set, "^[A-Z]+")]),
    shade = SHADES[2 - as.integer(block) %% 2]
  )

grid <- pathway_grid(
  transmute(literature,
    term = set, column = contrast, value = NES, p = list_fdr, selected = list_fdr <= ALPHA
  ),
  rows,
  tibble(
    column = COMPARISON_CONTRASTS, column_block = "Contrast",
    label = CONTRAST_SHORT[COMPARISON_CONTRASTS]
  ),
  title = "Literature-informed and ROS sets",
  subtitle = "Grouped by hypothesis; [expected direction]",
  tones = set_names(c("Contrast", levels(rows$block))),
  x_angle = 90, selected_label = paste("FDR ≤", ALPHA), fill_limit = 2.5, outline_shape = 21,
  min_block_mm = MIN_BLOCK_MM
) +
  theme(strip.text.y.left = element_text(
    angle = 0, hjust = 0, size = 5, margin = margin(1.5, 1.5, 1.5, 1.5)
  ))

height_mm <- sum(pmax(3 * table(rows$block), MIN_BLOCK_MM)) + nlevels(rows$block) + 30
save_figure(grid, "F03_Pathways", "S5b_literature", 116, height_mm)
list(
  plot = grid,
  data = literature |>
    left_join(select(rows, set = term, label, collection), by = "set") |>
    transmute(hypothesis, expected, set, label, collection, contrast, NES, p, fdr = list_fdr, size),
  height_mm = height_mm
)
