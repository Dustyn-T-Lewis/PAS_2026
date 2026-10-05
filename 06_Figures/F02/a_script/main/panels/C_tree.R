#!/usr/bin/env Rscript
# F02 c. The 15 strongest female and male terms, clustered on shared genes (08_Display), with their
# NES and FDR in each sex and for the treatment-by-sex test.

source(here::here("R", "panels.R"))

COLUMNS <- c("PAS_vs_VEH_female", "PAS_vs_VEH_male", "treatment_by_sex")

inputs <- c(display = here("03_Pathway_Enrichment", "08_Display", "c_data", "display.rds"))
display <- readRDS(inputs[["display"]])
clusters <- display$clusters |>
  mutate(cluster = if_else(
    size > 1, enrichVolcano::clean_label(cluster_name, width = 200), "Single term"
  ))
names <- sort(unique(clusters$cluster[clusters$size > 1]))
cluster_colours <- c(
  set_names(RColorBrewer::brewer.pal(length(names), "Dark2"), names),
  "Single term" = "grey85"
)
data <- display$tree_dots |>
  filter(contrast %in% COLUMNS) |>
  left_join(select(clusters, term, group = cluster), by = "term") |>
  mutate(
    column = factor(CONTRAST_LABELS[contrast], CONTRAST_LABELS[COLUMNS]),
    label = str_trunc(enrichVolcano::clean_label(term, width = 200), 46),
    padj = padj
  )

panel <- pathway_matrix(data,
  title = "Strongest female and male terms, clustered on shared genes",
  subtitle = paste(
    "F, M: PAS in each sex; F×M: interaction; clusters named by their top term;",
    "\nblack ring: FDR ≤ 0.05"
  ),
  tree = display$tree, group_colours = cluster_colours, group_name = "Cluster",
  column_names = set_names(c("F", "M", "F×M"), CONTRAST_LABELS[COLUMNS])
)
panel_data <- select(data, term, group, contrast, NES, p, padj)
panel_note <- "Tree terms with cluster, NES, p and FDR per contrast (08_Display)."
save_panel(panel, "F02", "main", "C_tree", 160, 170)
