#!/usr/bin/env Rscript
# Figure 3, whether the sexes respond to PAS differently: (a) female against male effect per
# protein, (b) sex-pattern labels per level, (c) signature carry-over between sexes, (d) interaction
# rings. Runs the panel scripts in panels/ and assembles them.

source(here::here("R", "panels.R"))

panels <- run_panels(
  "F03", "main", c("A_female_male", "B_labels", "C_signature", "D_interaction_rings")
)
plots <- map(panels, "plot")
tag <- \(p, letter) p + labs(tag = letter)
top <- tag(plots$A_female_male, "a") + tag(plots$B_labels, "b") + plot_layout(widths = c(1.2, 1))
left <- top / tag(plots$C_signature, "c") + plot_layout(heights = c(1.2, 1))
right <- wrap_elements(full = plots$D_interaction_rings) + labs(tag = "d")
figure <- (left | right) + plot_layout(widths = c(1.35, 1)) &
  theme(plot.tag = element_text(size = 9, face = "bold"))
save_figure(figure, figure_dirs("F03", "main")$reports, "F03",
  height_mm = 190, width_mm = LANDSCAPE_MM
)
write_panel_data(panels, "F03", "F03")

sessionInfo()
