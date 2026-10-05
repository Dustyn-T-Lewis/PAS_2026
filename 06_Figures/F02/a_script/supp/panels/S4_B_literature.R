#!/usr/bin/env Rscript
# S4 b. Literature-informed pathways (Burke 2026, Wang 2019) and the ROS sets in the three PAS
# contrasts, BH across the list (08_Display). Chosen after results; all shown.

source(here::here("R", "panels.R"))

inputs <- c(display = here("03_Pathway_Enrichment", "08_Display", "c_data", "display.rds"))
data <- readRDS(inputs[["display"]])$literature |>
  mutate(
    column = factor(CONTRAST_LABELS[contrast], CONTRAST_LABELS[PRIMARY]),
    label = str_glue("{str_trunc(enrichVolcano::clean_label(term, width = 200), 38)} ({expected})"),
    group = factor(hypothesis, unique(hypothesis)),
    padj = list_fdr
  ) |>
  arrange(group, term) |>
  mutate(term = factor(term, unique(term)))
hypothesis_colours <- set_names(
  RColorBrewer::brewer.pal(nlevels(data$group), "Set2"), levels(data$group)
)
panel <- pathway_matrix(data,
  title = "Literature-informed pathways",
  subtitle = paste(
    "F, M: PAS in each sex; F×M: interaction; (expected direction);",
    "\nsize and ring: FDR across the list"
  ),
  group_colours = hypothesis_colours, group_name = "Hypothesis",
  column_names = set_names(c("F", "M", "F×M"), CONTRAST_LABELS[PRIMARY])
)
panel_data <- select(data, contrast, hypothesis, expected, database, term, NES, p, padj = list_fdr)
panel_note <- "Literature and ROS sets: fgsea NES, p and BH across the list (08_Display)."
save_panel(panel, "F02", "supp", "S4_B_literature", 140, 110)
