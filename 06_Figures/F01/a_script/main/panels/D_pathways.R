#!/usr/bin/env Rscript
# F01 d. Pathways per contrast at fgsea FDR ≤ 0.05, every tested term: an up bar and a down bar
# per contrast, each stacked by database (a term sits in one database, so nothing counts twice).

source(here::here("R", "panels.R"))

SHADES <- list(
  Up = c("#67000D", "#A50F15", "#EF3B2C", "#FC9272", "#FCBBA1"),
  Down = c("#08306B", "#08519C", "#4292C6", "#9ECAE1", "#C6DBEF")
)
OFFSET <- c(Up = -0.2, Down = 0.2)

inputs <- c(enrichment = here("03_Pathway_Enrichment", "02_Contrasts", "c_data", "enrichment.rds"))
databases <- names(DATABASE_COLOURS)
counts <- readRDS(inputs[["enrichment"]]) |>
  map(\(tests) tests$fgsea@results) |>
  list_rbind() |>
  filter(padj <= ALPHA) |>
  count(contrast, database, direction = if_else(score > 0, "Up", "Down")) |>
  mutate(
    contrast = factor(CONTRAST_LABELS[contrast], CONTRAST_LABELS),
    x = as.numeric(contrast) + OFFSET[direction],
    database = factor(database, databases),
    fill = factor(
      str_glue("{direction}: {database}"),
      c(str_glue("Up: {databases}"), str_glue("Down: {databases}"))
    )
  )
totals <- summarise(counts, n = sum(n), .by = c(contrast, direction, x))
bands <- tibble(contrast = factor(CONTRAST_LABELS, CONTRAST_LABELS)) |>
  mutate(x = as.numeric(contrast), band = CONTRAST_COLOURS[as.character(contrast)])

panel <- ggplot(counts, aes(x, n)) +
  geom_tile(
    data = bands, aes(x = x, y = 0, fill = NULL), fill = bands$band, width = 1, height = Inf,
    alpha = 0.12, inherit.aes = FALSE
  ) +
  geom_col(aes(fill = fill), width = 0.38, colour = "white", linewidth = 0.15) +
  geom_text(data = totals, aes(label = n), vjust = -0.4, size = 1.9) +
  scale_fill_manual(
    values = set_names(c(SHADES$Up, SHADES$Down), levels(counts$fill)),
    name = NULL, drop = FALSE
  ) +
  scale_x_continuous(breaks = seq_along(CONTRAST_LABELS), labels = str_replace(
    CONTRAST_LABELS, ", ", ",\n"
  ), expand = expansion(add = 0.1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = NULL, y = "Pathways at FDR ≤ 0.05", title = "Pathways per contrast",
    subtitle = "fgsea, every tested term; left bar up (reds), right bar down (blues), by database"
  ) +
  theme_figure() +
  theme(
    panel.grid.major.x = element_blank(), legend.key.size = unit(2.3, "mm"),
    legend.text = element_text(size = 5)
  ) +
  guides(fill = guide_legend(ncol = 2))
panel_data <- select(counts, contrast, direction, database, terms = n)
panel_note <- "fgsea terms at FDR ≤ 0.05 per contrast, direction and database, every tested term."
save_panel(panel, "F01", "main", "D_pathways", 120, 80)
