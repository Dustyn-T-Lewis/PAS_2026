#!/usr/bin/env Rscript
# F04 b. Proteins and pathways that track each phenotype: pooled, female-only and slope difference.

source(here::here("R", "panels.R"))

PHENOTYPE_LABELS <- c(
  csa_mean = "CSA", capillary_density = "Capillary density",
  satellite_cells_per_fiber = "Satellite cells / fibre",
  conductance_sol_pyr = "Soleus pyruvate conductance",
  conductance_pla_pyr = "Plantaris pyruvate conductance",
  conductance_sol_oct = "Soleus octanoyl conductance"
)
MODEL_LABELS <- c(
  "pooled slope" = "Pooled", "females only" = "Female",
  "slope difference, male minus female" = "Δ slope"
)
INPUT_LABELS <- c(
  "phenotypes only" = "Phenotypes", "proteins features" = "Proteins",
  "sets features" = "Pathways", "modules features" = "Modules",
  "combined features" = "Proteome combined",
  "combined features + phenotypes" = "Proteome + phenotypes"
)
SEX_COLOURS_SHORT <- c(F = "#D55E00", M = "#0072B2")

inputs <- c(
  proteins = here("02_Differential_Expression", "04_Associate", "c_data", "04_associate.xlsx"),
  pathways = here("03_Pathway_Enrichment", "04_Associate", "c_data", "04_associate.xlsx")
)
phenotype_factor <- \(x) factor(PHENOTYPE_LABELS[x], rev(PHENOTYPE_LABELS))

protein_counts <- readxl::read_excel(inputs[["proteins"]], "summary") |>
  transmute(phenotype, model, level = "Proteins\n(FDR ≤ 0.10)", n = fdr_le_0.10)
pathway_counts <- readxl::read_excel(inputs[["pathways"]], "phenotype_gsea", guess_max = 1e5) |>
  summarise(n = sum(padj <= ALPHA), .by = c(phenotype, model)) |>
  mutate(level = "Pathways\n(fgsea FDR ≤ 0.05)")
counts <- bind_rows(protein_counts, pathway_counts) |>
  mutate(
    phenotype = phenotype_factor(phenotype),
    model = factor(MODEL_LABELS[model], MODEL_LABELS),
    level = factor(level, c("Proteins\n(FDR ≤ 0.10)", "Pathways\n(fgsea FDR ≤ 0.05)"))
  )
panel <- ggplot(counts, aes(model, phenotype, fill = log1p(n))) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = n, colour = log1p(n) > 3.5), size = 2) +
  facet_wrap(~level) +
  scale_fill_distiller(palette = "Greys", direction = 1, guide = "none") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none") +
  labs(
    x = NULL, y = NULL, title = "How much of the proteome tracks each phenotype",
    subtitle = "Pooled: within-group slope, 20 mice; female: 10 mice; Δ slope: male minus female"
  ) +
  theme_figure() +
  theme(axis.ticks = element_blank(), panel.grid = element_blank())

panel_data <- select(counts, phenotype, model, level, n)
panel_note <- "Proteins (FDR ≤ 0.10) and pathways (FDR ≤ 0.05) per phenotype and model."
save_panel(panel, "F04", "main", "B_proteome_links", 130, 75)
