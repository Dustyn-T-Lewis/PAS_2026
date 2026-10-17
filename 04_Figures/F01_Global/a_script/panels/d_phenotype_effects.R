# F01 d: what PAS did to each phenotype within each sex, as a % change with its 95% CI on a log
# ratio axis (so rises and falls of the same size look the same size), with the treatment-by-sex p
# per trait.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(fit = output_path("02_Statistics", "01_Fit", "fit.rds"))
fit <- readRDS(inputs[["fit"]])
interaction_p <- with(fit$phenotype_interaction, set_names(P.Value, PHENOTYPE_LABELS[phenotype]))
effects <- fit$phenotype_effects |>
  mutate(
    sex = factor(if_else(sex == "F", "Female", "Male"), c("Male", "Female")),
    trait = factor(PHENOTYPE_LABELS[phenotype], rev(PHENOTYPE_LABELS)),
    across(c(pct_change, ci_low, ci_high), \(x) 1 + x / 100)
  )

plot <- ggplot(effects, aes(pct_change, trait, colour = sex)) +
  geom_vline(xintercept = 1, colour = "grey50", linewidth = 0.3) +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high, shape = P.Value <= ALPHA),
    position = position_dodge(width = 0.6), size = 0.25, linewidth = 0.4, fill = "white"
  ) +
  scale_x_continuous(
    transform = "log2", breaks = c(0.5, 1, 2, 4), labels = \(x) sprintf("%+.0f%%", 100 * (x - 1))
  ) +
  scale_y_discrete(labels = \(x) {
    sprintf("%s\n%s p %.2g", x, CONTRAST_SHORT[["treatment_by_sex"]], interaction_p[x])
  }) +
  scale_colour_manual(values = SEX_COLOURS, breaks = c("Female", "Male"), name = NULL) +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 21),
    labels = c(`TRUE` = paste("filled: p ≤", ALPHA), `FALSE` = paste("open: p >", ALPHA)),
    name = NULL
  ) +
  labs(
    x = "% change (95% CI, log scale)", y = NULL, title = "Phenotype effects",
    subtitle = "PAS vs vehicle within sex"
  ) +
  theme_figure() +
  theme(legend.box = "horizontal", legend.spacing.x = unit(2, "mm")) +
  guides(
    colour = guide_legend(order = 1, ncol = 1),
    shape = guide_legend(order = 2, ncol = 1, reverse = TRUE)
  )

save_panel(plot, "F01_Global", "d_phenotype_effects", 66, 64)
list(
  plot = plot,
  data = left_join(fit$phenotype_effects,
    select(fit$phenotype_interaction, phenotype, interaction_p = P.Value),
    by = "phenotype"
  ),
  inputs = inputs
)
