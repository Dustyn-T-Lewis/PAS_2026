#!/usr/bin/env Rscript
# F03 b. Features at p ≤ 0.05 per sex-pattern label and level (05_Summary/04_Sex_Patterns).

source(here::here("R", "panels.R"))

LABEL_COLOURS <- c(
  "sex-differential" = "#7B3294", shared = "#009E73", "detected in females" = "#D55E00",
  "detected in males" = "#0072B2", neither = "grey85"
)
SEX_NAMES <- c(PAS_vs_VEH_female = "Female", PAS_vs_VEH_male = "Male")
COLLECTIONS <- list(full = c("Hallmark", "Reactome", "KEGG", "GO:BP"), mito = "MitoCarta")
UNIVERSE_LABELS <- c(full = "full proteome", mito = "mitochondrial proteome")
GROUP_ORDER <- c("Proteins, full", "Proteins, mito", "Pathways, full", "Pathways, mito", "Modules")
N_LABELLED <- 8

inputs <- c(
  patterns = here("05_Summary", "04_Sex_Patterns", "c_data", "04_sex_patterns.xlsx")
)
patterns <- readxl::read_excel(inputs[["patterns"]], "patterns", guess_max = 1e5) |>
  mutate(label = factor(label, names(LABEL_COLOURS)))

counts <- patterns |>
  filter(label != "neither") |>
  mutate(group = case_when(
    level == "protein" ~ str_glue("Proteins, {universe}"),
    level == "pathway" ~ str_glue("Pathways, {universe}"),
    TRUE ~ "Modules"
  )) |>
  count(group, label) |>
  complete(group = GROUP_ORDER, label, fill = list(n = 0)) |>
  mutate(group = factor(group, rev(GROUP_ORDER)))
totals <- summarise(counts, n = sum(n), .by = group)
panel <- ggplot(counts, aes(n, group)) +
  geom_col(aes(fill = label), width = 0.7) +
  geom_text(data = totals, aes(label = n), hjust = -0.2, size = 1.9) +
  scale_fill_manual(values = LABEL_COLOURS, breaks = names(LABEL_COLOURS)[1:4], name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "Features at p ≤ 0.05", y = NULL, title = "How the response divides by sex",
    subtitle = "One label per feature, so bars stack; no module reaches p ≤ 0.05"
  ) +
  theme_figure() +
  theme(legend.position = "inside", legend.position.inside = c(0.75, 0.25))

panel_data <- counts
panel_note <- "Features per label and level at p ≤ 0.05 (05_Summary/04_Sex_Patterns)."
save_panel(panel, "F03", "main", "B_labels", 85, 85)
