#!/usr/bin/env Rscript
# S2 b. Distance of each mouse from its sex's vehicle centroid; RRPP test.

source(here::here("R", "panels.R"))

inputs <- c(magnitude = here("05_Summary", "06_Magnitude", "c_data", "magnitude.rds"))
magnitude <- readRDS(inputs[["magnitude"]])
stat <- \(t) with(filter(magnitude$tests, term == t), sprintf("R² %.2f, p = %.2f", r2, p))
panel <- magnitude$magnitude |>
  mutate(group = factor(group, FIGURE_GROUPS)) |>
  ggplot(aes(group, distance, colour = group, shape = sex)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.25, colour = "grey40") +
  geom_point(position = position_jitter(width = 0.1, seed = SEED), size = 1.5) +
  scale_colour_manual(values = GROUP_COLOURS, guide = "none") +
  scale_shape_manual(values = SEX_SHAPES, guide = "none") +
  labs(
    x = NULL, y = "Distance from own-sex vehicle centroid", title = "Response magnitude",
    subtitle = str_glue(
      "RRPP: treatment {stat('treatment')}; sex × treatment {stat('sex:treatment')};\n",
      "PAS moves neither sex's proteome as a whole beyond vehicle spread"
    )
  ) +
  theme_figure()
panel_data <- magnitude$tests
panel_note <- "RRPP ANOVA on distance from own-sex vehicle centroid (05_Summary/06_Magnitude)."
save_panel(panel, "F01", "supp", "S2_B_magnitude", 80, 65)
