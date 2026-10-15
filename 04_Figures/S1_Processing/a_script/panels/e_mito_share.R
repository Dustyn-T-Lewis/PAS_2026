# S1 e: the share of each pellet's protein signal from MitoCarta proteins, by group, with the
# sex and treatment p from 02_Quantification.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(
  quantification = output_path("01_Processing", "02_Quantification", "02_quantification.xlsx")
)
share <- readxl::read_excel(inputs[["quantification"]], "mito_share") |>
  mutate(group = factor(group, FIGURE_GROUPS))
tests <- readxl::read_excel(inputs[["quantification"]], "mito_share_tests")
label <- paste(
  sprintf("%s p = %.2g", str_to_sentence(tests$factor), tests$p_value),
  collapse = "\n"
)

plot <- ggplot(share, aes(group, pct_mito_signal, colour = group)) +
  geom_point(position = position_jitter(width = 0.15, seed = SEED), size = 1.3) +
  annotate("label",
    x = -Inf, y = Inf, label = label, hjust = 0, vjust = 1, size = 1.7, lineheight = 0.95,
    fill = NA, linewidth = 0, label.padding = unit(1, "mm")
  ) +
  scale_colour_manual(values = GROUP_COLOURS, guide = "none") +
  scale_x_discrete(labels = GROUP_SHORT) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  labs(
    x = NULL, y = "% signal from MitoCarta", title = "Pellet mitochondrial share",
    subtitle = "Differs by sex, not by treatment"
  ) +
  theme_figure()

save_panel(plot, "S1_Processing", "e_mito_share", 60, 62)
list(plot = plot, data = share, inputs = inputs)
