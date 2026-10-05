#!/usr/bin/env Rscript
# F04 c. Female PAS signature score against the phenotypes PAS changed in females.

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

changed <- names(PHENOTYPE_LABELS)[1:3]
slopes <- sp$slopes |>
  filter(signature == "female", phenotype %in% changed) |>
  select(phenotype, model, p) |>
  pivot_wider(names_from = model, values_from = p)
facet_labels <- set_names(
  str_glue(
    "{PHENOTYPE_LABELS[slopes$phenotype]}\nslope p = {round(slopes$`pooled slope`, 2)}; ",
    "Δ p = {round(slopes$`slope difference, male minus female`, 2)}"
  ),
  slopes$phenotype
)
ph <- read_phenotypes(sp$scores$sample_id[sp$scores$signature == "female"])
score_data <- sp$scores |>
  filter(signature == "female") |>
  bind_cols(select(ph, all_of(changed))) |>
  pivot_longer(all_of(changed), names_to = "phenotype", values_to = "value") |>
  mutate(phenotype = factor(facet_labels[phenotype], facet_labels))
panel <- ggplot(score_data, aes(value, score, colour = group, shape = sex)) +
  geom_point(size = 1.7) +
  facet_wrap(~phenotype, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = GROUP_COLOURS, name = NULL) +
  scale_shape_manual(values = SEX_SHAPES, guide = "none") +
  labs(
    x = "Phenotype value", y = "Female PAS signature score",
    title = "The female PAS signature against the phenotypes PAS changed in females",
    subtitle = paste(
      "singscore on the female signature proteins; slope: within-group (pooled); Δ: male minus",
      "female.\nGroups separate on both axes, but within groups the score does not track the trait"
    )
  ) +
  theme_figure() +
  theme(strip.text = element_text(size = 5.5, face = "plain"))

panel_data <- filter(sp$slopes, signature == "female")
panel_note <- "Female signature score against each phenotype: pooled, female, slope difference."
save_panel(panel, "F04", "main", "C_signature_score", 150, 75)
