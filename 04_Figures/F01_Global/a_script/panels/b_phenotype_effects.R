# F01 b: what PAS did to each phenotype within each sex, with the treatment-by-sex p per trait.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(fit = output_path("02_Statistics", "01_Fit", "fit.rds"))
fit <- readRDS(inputs[["fit"]])
interaction_p <- with(fit$phenotype_interaction, set_names(P.Value, PHENOTYPE_LABELS[phenotype]))

plot <- phenotype_effect_plot(fit$phenotype_effects,
  title = "Phenotype effects",
  subtitle = "PAS vs vehicle within sex; F − M: interaction p"
) +
  scale_y_discrete(labels = \(x) sprintf("%s\nF − M p %.2g", x, interaction_p[x])) +
  scale_x_continuous(
    transform = "log2", breaks = c(0.5, 1, 2, 4), labels = \(x) sprintf("%+.0f%%", 100 * (x - 1))
  ) +
  labs(x = "% change (95% CI, log scale)") +
  theme(legend.position = "bottom", legend.key.size = unit(2.5, "mm"), legend.box = "vertical") +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))

save_panel(plot, "F01_Global", "b_phenotype_effects", 54, 62)
list(
  plot = plot,
  data = left_join(fit$phenotype_effects,
    select(fit$phenotype_interaction, phenotype, interaction_p = P.Value),
    by = "phenotype"
  ),
  inputs = inputs
)
