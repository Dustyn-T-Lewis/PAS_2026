#!/usr/bin/env Rscript
# F03 c. Whether each sex's PAS signature moves the same way in the other sex (limma::fry,
# 05_Summary/07_Signature).

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
  signature = here("05_Summary", "07_Signature", "c_data", "signature.rds")
)
signature <- readRDS(inputs[["signature"]])

panel_name <- \(source, target) {
  str_glue("{SEX_NAMES[source]} signature in the {tolower(SEX_NAMES[target])} ranking")
}
rank_data <- signature$ranks |>
  mutate(
    panel = panel_name(source, target),
    direction = factor(str_to_sentence(direction), c("Up", "Down"))
  )
rank_tests <- signature$tests |>
  mutate(
    panel = panel_name(source, target),
    direction = factor(str_to_sentence(direction), c("Up", "Down")),
    label = sprintf(
      "%s: %d proteins, %.0f%% same way, fry p = %.2f", direction, proteins, 100 * same_way, p
    )
  )
panel <- ggplot(rank_data, aes(rank, direction, colour = direction)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3) +
  geom_point(shape = "|", size = 3) +
  geom_label(
    data = rank_tests, aes(x = 0, label = label), hjust = 0, vjust = -0.7, size = 1.8,
    colour = "grey20", fill = "white", linewidth = 0, label.padding = unit(0.6, "mm")
  ) +
  facet_wrap(~panel, ncol = 1) +
  scale_colour_manual(values = c(Up = "#D55E00", Down = "#2166AC"), guide = "none") +
  scale_x_continuous(
    limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("most up", "", "most down")
  ) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.2))) +
  labs(
    x = "Rank in the other sex's PAS contrast", y = NULL,
    title = "Does one sex's signature carry over?",
    subtitle = "Proteins at p ≤ 0.05 for PAS in one sex, placed in the other sex's ranking"
  ) +
  theme_figure() +
  theme(strip.text = element_text(size = 6))

panel_data <- signature$tests
panel_note <- "fry test of each sex's PAS signature in the other sex (05_Summary/07_Signature)."
save_panel(panel, "F03", "main", "C_signature", 180, 75)
