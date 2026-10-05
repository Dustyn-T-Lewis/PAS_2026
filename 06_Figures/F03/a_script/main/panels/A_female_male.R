#!/usr/bin/env Rscript
# F03 a. Female against male PAS effect per protein, coloured by the interaction-anchored sex
# pattern (05_Summary/04_Sex_Patterns).

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

proteins <- filter(patterns, level == "protein", universe == "full")
panel <- ggplot(arrange(proteins, label != "neither"), aes(effect_m, effect_f, colour = label)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_abline(linetype = "dotted", linewidth = 0.3) +
  geom_point(size = 0.6) +
  ggrepel::geom_text_repel(
    data = slice_min(proteins, p_int, n = N_LABELLED), aes(label = gene), size = 1.9,
    colour = "black", min.segment.length = 0, segment.size = 0.2, max.overlaps = Inf,
    seed = SEED
  ) +
  scale_colour_manual(values = LABEL_COLOURS, guide = "none") +
  labs(
    x = "log2 FC, PAS vs VEH, males", y = "log2 FC, PAS vs VEH, females",
    title = "Female against male PAS effect",
    subtitle = str_glue(
      "Full proteome, {nrow(proteins)} proteins; labels from p ≤ 0.05 in each sex and the ",
      "interaction;\nlabelled: {N_LABELLED} smallest interaction p"
    )
  ) +
  theme_figure()

panel_data <- select(
  proteins, feature, gene, label, effect_f, effect_m, effect_int, p_f, p_m, p_int, fdr_int
)
panel_note <- "Female and male PAS effect per protein, full universe, with the sex-pattern label."
save_panel(panel, "F03", "main", "A_female_male", 95, 85)
