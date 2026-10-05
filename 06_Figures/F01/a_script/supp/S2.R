#!/usr/bin/env Rscript
# Supplementary Figure 2, further global views: variance by component, response magnitude, effect
# sizes, the mitochondrial universe, female against male effects, p-value distributions and the
# baseline sex pathways. Runs the panel scripts in panels/ and assembles them.

source(here::here("R", "panels.R"))

panels <- run_panels("F01", "supp", c(
  "S2_A_variance", "S2_B_magnitude", "S2_C_effect_sizes", "S2_D_mito_mds", "S2_E_female_male",
  "S2_F_mito_p", "S2_G_sex_p", "S2_H_sex_pathways"
))
plots <- map(panels, "plot")
row_1 <- plots$S2_A_variance + plots$S2_B_magnitude + plots$S2_C_effect_sizes +
  plot_layout(widths = c(1.1, 0.9, 1))
row_2 <- plots$S2_D_mito_mds + plots$S2_E_female_male + plot_layout(widths = c(1, 1.4))
row_3 <- plots$S2_F_mito_p + plots$S2_G_sex_p + plots$S2_H_sex_pathways +
  plot_layout(widths = c(1.3, 0.9, 1.1))
figure <- row_1 / row_2 / row_3 + plot_annotation(tag_levels = "a")
save_figure(figure, figure_dirs("F01", "supp")$reports, "S2",
  height_mm = 230, width_mm = LANDSCAPE_MM
)
write_panel_data(panels, "F01", "S2")

sessionInfo()
