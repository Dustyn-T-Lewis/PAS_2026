#!/usr/bin/env Rscript
# S2 a. R² of sex, treatment, their interaction and the MitoCarta share on PC1 to PC4.

source(here::here("R", "panels.R"))

inputs <- c(global = here("05_Summary", "05_Global", "c_data", "global.rds"))
component_r2 <- readRDS(inputs[["global"]])$component_r2
factors <- unique(component_r2$factor)
panel <- ggplot(component_r2, aes(pc, factor(factor, rev(factors)), fill = r2)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", r2), colour = r2 > 0.5), size = 1.9) +
  scale_fill_distiller(palette = "Greys", direction = 1, limits = c(0, 1), name = "R²") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"), guide = "none") +
  labs(
    x = NULL, y = NULL, title = "Variance by component",
    subtitle = "R² of each factor alone; PC1 is sex and MitoCarta share together"
  ) +
  theme_figure() +
  theme(axis.ticks = element_blank(), panel.grid = element_blank())
panel_data <- select(component_r2, component, factor, r2)
panel_note <- "R² of each factor and the MitoCarta share on PC1 to PC4 (05_Summary/05_Global)."
save_panel(panel, "F01", "supp", "S2_A_variance", 90, 60)
