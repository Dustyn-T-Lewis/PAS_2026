# S1 d: the surrogate variable (sva, one SV) every model carries, against each run's depth, with
# its r per sex and the treatment p from 03_Design.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(design = output_path("01_Processing", "03_Design", "03_design.xlsx"))
sv <- readxl::read_excel(inputs[["design"]], "surrogate") |>
  mutate(group = factor(group, FIGURE_GROUPS))
r <- readxl::read_excel(inputs[["design"]], "surrogate_r")
label <- with(r, paste(
  sprintf(
    "%s r = %.2f%s", str_to_sentence(mice), r_with_run_depth,
    if_else(is.na(treatment_p), "", sprintf(", treatment p = %.2f", treatment_p))
  ),
  collapse = "\n"
))

plot <- ggplot(sv, aes(run_depth_log2, sv1)) +
  geom_point(aes(colour = group, shape = sex), size = 1.3) +
  annotate("label",
    x = -Inf, y = Inf, label = label, hjust = 0, vjust = 1, size = 1.7, lineheight = 0.95,
    fill = NA, linewidth = 0, label.padding = unit(1, "mm")
  ) +
  scale_colour_manual(
    values = GROUP_COLOURS, labels = GROUP_LABELS, name = NULL,
    guide = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  scale_shape_manual(values = SEX_SHAPES, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.45))) +
  labs(
    x = "run depth (log2 precursor signal)", y = "surrogate variable 1",
    title = "Surrogate variable", subtitle = "Tracks run depth, not treatment"
  ) +
  theme_figure()

save_panel(plot, "S1_Processing", "d_surrogate", 60, 62)
list(plot = plot, data = sv, inputs = inputs)
