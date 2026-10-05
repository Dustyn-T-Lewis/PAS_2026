#!/usr/bin/env Rscript
# S2 h. Strongest baseline sex differences in pathways, male minus female (08_Display).

source(here::here("R", "panels.R"))

SEX_COLUMNS <- c(sex_in_VEH = "Vehicle", sex_in_PAS = "PAS")
inputs <- c(display = here("03_Pathway_Enrichment", "08_Display", "c_data", "display.rds"))
sex_data <- readRDS(inputs[["display"]])$sex_dots |>
  mutate(
    column = factor(SEX_COLUMNS[contrast], SEX_COLUMNS),
    label = str_trunc(enrichVolcano::clean_label(term, width = 200), 40),
    label = reorder(label, if_else(contrast == "sex_in_VEH", NES, 0))
  )
panel <- ggplot(sex_data, aes(column, label)) +
  geom_point(aes(size = -log10(padj), fill = NES), shape = 21, stroke = 0.25) +
  scale_fill_distiller(palette = "PuOr", direction = -1, limits = c(-2.6, 2.6)) +
  scale_size_area(max_size = 3.5, name = "-log10 FDR") +
  scale_x_discrete(expand = expansion(add = 0.7)) +
  labs(
    x = NULL, y = NULL, fill = "NES (> 0:\nhigher in males)",
    title = "Baseline sex differences in pathways",
    subtitle = "10 strongest collapsed terms per contrast"
  ) +
  theme_figure() +
  theme(axis.text.y = element_text(size = 5))
panel_data <- select(sex_data, contrast, database, term, NES, p, padj)
panel_note <- "Baseline sex terms, male minus female, in vehicle and in PAS (08_Display)."
save_panel(panel, "F01", "supp", "S2_H_sex_pathways", 100, 90)
