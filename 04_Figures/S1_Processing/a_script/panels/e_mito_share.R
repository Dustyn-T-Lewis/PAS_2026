# S1 e: the share of each pellet's protein signal from MitoCarta proteins, by group, with the sex
# p from 02_Quantification and the PAS effect within each sex (difference, 95% CI and Welch p)
# from 05_Composition.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(
  quantification = output_path("01_Processing", "02_Quantification", "proteins.rds"),
  composition = output_path("02_Statistics", "05_Composition", "composition.rds")
)
quant <- readRDS(inputs[["quantification"]])
share <- quant$mito_share |>
  mutate(group = factor(group, FIGURE_GROUPS))
sex_p <- quant$mito_share_tests |>
  filter(factor == "sex") |>
  pull(p_value)
composition <- readRDS(inputs[["composition"]])
within <- filter(composition$tests, group == composition$groups$mito_share)
label <- paste(
  c(
    sprintf("Sex p = %.2g", sex_p),
    with(within, sprintf(
      "PAS in %ss %+.1f points (%.1f to %.1f), p = %.2f", sex, difference, ci_low, ci_high, p
    ))
  ),
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
    subtitle = "One point per mouse"
  ) +
  theme_figure()

save_panel(plot, "S1_Processing", "e_mito_share", 60, 62)
list(plot = plot, data = share, inputs = inputs)
