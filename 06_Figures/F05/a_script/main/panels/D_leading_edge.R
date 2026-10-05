#!/usr/bin/env Rscript
# F05 d. Leading-edge proteins of the strongest interaction pathways, PAS log2 fold change per sex.

source(here::here("R", "panels.R"))

PHENOTYPE_SHORT <- c(
  conductance_sol_pyr = "Soleus pyruvate", capillary_density = "Capillary density",
  csa_mean = "CSA", conductance_sol_oct = "Soleus octanoyl"
)
MODEL_SHORT <- c("females only" = "F", "slope difference, male minus female" = "Δ")
term_label <- \(x, n = 40) str_trunc(enrichVolcano::clean_label(x, width = 200), n)
N_EDGE_TERMS <- 6
N_EDGE_GENES <- 5

inputs <- c(
  display = here("03_Pathway_Enrichment", "08_Display", "c_data", "display.rds"),
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)

display <- readRDS(inputs[["display"]])
de <- readRDS(inputs[["contrasts"]])
edge_terms <- display$ring_terms |>
  filter(universe == "full", contrast == "treatment_by_sex") |>
  slice_min(p, n = N_EDGE_TERMS, with_ties = FALSE) |>
  mutate(
    gene = map(leading_edge, \(x) head(unlist(x), N_EDGE_GENES)), label = term_label(term, 30)
  ) |>
  select(term, label, gene) |>
  unnest(gene)
edge_effects <- de$results |>
  filter(
    universe == "full", contrast %in% c("PAS_vs_VEH_female", "PAS_vs_VEH_male"),
    gene %in% edge_terms$gene
  ) |>
  slice_max(AveExpr, n = 1, by = c(gene, contrast), with_ties = FALSE) |>
  transmute(gene,
    column = if_else(contrast == "PAS_vs_VEH_female", "Female", "Male"), logFC,
    p = P.Value
  ) |>
  inner_join(edge_terms, by = "gene", relationship = "many-to-many") |>
  mutate(label = factor(label, unique(edge_terms$label)))
panel <- ggplot(edge_effects, aes(column, gene, fill = logFC)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(data = \(d) filter(d, p <= ALPHA), label = "•", size = 2.2) +
  facet_grid(label ~ ., scales = "free_y", space = "free_y") +
  scale_fill_distiller(
    palette = "RdBu", limits = c(-1.2, 1.2), oob = scales::squish,
    name = "log2 FC"
  ) +
  labs(
    x = NULL, y = NULL, title = "Proteins behind the sex-different pathways",
    subtitle = str_glue(
      "Leading edge (first {N_EDGE_GENES}) of the {N_EDGE_TERMS} strongest interaction terms;\n",
      "PAS log2 FC per sex; • p ≤ 0.05"
    )
  ) +
  theme_figure() +
  theme(
    axis.text.y = element_text(size = 5),
    strip.text.y = element_text(size = 5, angle = 0, hjust = 0),
    panel.grid = element_blank(), panel.spacing = unit(0.6, "mm")
  )

panel_data <- select(edge_effects, label, gene, column, logFC, p)
panel_note <- "PAS log2 FC per sex, leading edge of the strongest interaction terms."
save_panel(panel, "F05", "main", "D_leading_edge", 110, 110)
