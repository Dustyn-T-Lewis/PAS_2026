#!/usr/bin/env Rscript
# F04 a. PAS minus vehicle for each phenotype in each sex (05_Summary/08_Signature_Phenotype).

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
  signature_phenotype = here(
    "05_Summary", "08_Signature_Phenotype", "c_data", "signature_phenotype.rds"
  )
)
sp <- readRDS(inputs[["signature_phenotype"]])
phenotype_factor <- \(x) factor(PHENOTYPE_LABELS[x], rev(PHENOTYPE_LABELS))

effects <- sp$phenotype_effects |>
  mutate(
    trait = if_else(interaction_p <= ALPHA, paste(PHENOTYPE_LABELS[phenotype], "†"),
      PHENOTYPE_LABELS[phenotype]
    ),
    trait = factor(trait, rev(unique(trait[order(match(phenotype, names(PHENOTYPE_LABELS)))]))),
    sex = factor(sex, c("M", "F"))
  )
panel <- ggplot(effects, aes(difference, trait, colour = sex)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_errorbar(aes(xmin = low, xmax = high),
    width = 0, linewidth = 0.5, orientation = "y", position = position_dodge(width = 0.55)
  ) +
  geom_point(aes(shape = p <= ALPHA), size = 1.8, position = position_dodge(width = 0.55)) +
  scale_colour_manual(
    values = SEX_COLOURS_SHORT, labels = c(F = "Female", M = "Male"),
    breaks = c("F", "M"), name = NULL
  ) +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1), labels = c("p > 0.05", "p ≤ 0.05"),
    name = NULL
  ) +
  labs(
    x = "PAS − vehicle (SD units, 95% CI)", y = NULL, title = "PAS effect on each phenotype",
    subtitle = "Welch t within sex; † treatment × sex p ≤ 0.05 (lm)"
  ) +
  theme_figure()

panel_data <- sp$phenotype_effects
panel_note <- "PAS minus vehicle per phenotype and sex (05_Summary/08_Signature_Phenotype)."
save_panel(panel, "F04", "main", "A_phenotype_effects", 120, 75)
