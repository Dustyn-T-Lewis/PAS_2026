#!/usr/bin/env Rscript
# S3 a. Compartment shifts by fry, a self-contained test of whether each GO Slim compartment's
# protein set moves as a whole, for PAS in each sex and for the treatment-by-sex interaction.

source(here::here("R", "panels.R"))

inputs <- c(
  fry = here("03_Pathway_Enrichment", "07_Compartments", "c_data", "compartments_fry.rds")
)
fry <- readRDS(inputs[["fry"]])@results |>
  filter(contrast %in% PRIMARY) |>
  mutate(
    contrast = factor(CONTRAST_LABELS[contrast], CONTRAST_LABELS[PRIMARY]),
    signed = if_else(direction == "up", 1, -1) * -log10(p),
    compartment = str_to_sentence(str_replace_all(str_remove(term, "^GOSLIM_"), "_", " "))
  ) |>
  mutate(order = signed[contrast == CONTRAST_LABELS[["PAS_vs_VEH_female"]]], .by = compartment) |>
  mutate(compartment = reorder(compartment, order))
cut <- -log10(ALPHA)

panel <- ggplot(fry, aes(signed, compartment, fill = contrast)) +
  geom_vline(xintercept = c(-cut, cut), linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_col(width = 0.7) +
  geom_text(data = \(d) filter(d, padj <= ALPHA), aes(label = "*"), hjust = -0.3, size = 2.5) +
  facet_wrap(~contrast, nrow = 1) +
  scale_fill_manual(values = CONTRAST_COLOURS, guide = "none") +
  labs(
    x = "Signed −log10 p (fry; > 0 up)", y = NULL, title = "Compartment shifts",
    subtitle = str_glue(
      "GO Slim cellular components; dashed: p = 0.05; * FDR ≤ 0.05 (none reach it, lowest ",
      "{signif(min(fry$padj), 2)})"
    )
  ) +
  theme_figure() +
  theme(axis.text.y = element_text(size = 5.5))
panel_data <- select(fry, contrast, term, size, direction, p, padj)
panel_note <- "fry per PAS contrast and GO Slim compartment (07_Compartments)."
save_panel(panel, "F01", "supp", "S3_A_compartments", 180, 90)
