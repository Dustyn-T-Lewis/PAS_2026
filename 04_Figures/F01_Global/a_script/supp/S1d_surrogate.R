# S1 d: the surrogate variable (sva, one SV) every model carries, against each run's depth, with
# its r per sex and the treatment p from 03_Design.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(design = output_path("01_Processing", "03_Design", "design.rds"))
design <- readRDS(inputs[["design"]])
sv <- design$surrogate |>
  mutate(group = factor(group, FIGURE_GROUPS))
r <- design$surrogate_r
label <- with(r, paste(
  sprintf(
    "%s r = %.2f%s", str_to_sentence(mice), r_with_run_depth,
    if_else(is.na(treatment_p), "", sprintf(", treatment p = %.2f", treatment_p))
  ),
  collapse = "\n"
))

plot <- ggplot(sv, aes(run_depth_log2, sv1)) +
  group_points() +
  annotate("label",
    x = -Inf, y = Inf, label = label, hjust = 0, vjust = 1, size = 1.7, lineheight = 0.95,
    fill = NA, linewidth = 0, label.padding = unit(1, "mm")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.45))) +
  labs(
    x = "run depth (log2 precursors identified)", y = "surrogate variable 1",
    title = "Surrogate variable", subtitle = "One point per mouse"
  ) +
  theme_figure()

save_figure(plot, "F01_Global", "S1d_surrogate", 96, 62)
list(plot = plot, data = sv)
