# S2 a: how much of each of the first four proteome principal components sex, treatment, their
# interaction, the pellet's MitoCarta share and run depth explain, as R² (04_PCA).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(pca = output_path("02_Statistics", "04_PCA", "pca.rds"))
r2 <- readRDS(inputs[["pca"]])$component_r2 |>
  mutate(
    pc = factor(pc, unique(pc)),
    factor = factor(factor, rev(unique(factor)))
  )

plot <- ggplot(r2, aes(pc, factor)) +
  geom_tile(aes(fill = r2), colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", r2), colour = r2 > 0.5), size = 1.8) +
  scale_fill_gradient(
    low = "grey97", high = "grey25", limits = c(0, 1), breaks = c(0, 0.5, 1), name = "R²",
    guide = colour_bar()
  ) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_x_discrete(position = "top", expand = expansion(0)) +
  scale_y_discrete(expand = expansion(0)) +
  labs(
    x = NULL, y = NULL, title = "What the proteome PCs track",
    subtitle = "R² of each factor on each component"
  ) +
  theme_figure() +
  theme(panel.grid = element_blank())

save_figure(plot, "F01_Global", "S2a_component_r2", 90, 55)
list(plot = plot, data = select(r2, -component))
