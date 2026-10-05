#!/usr/bin/env Rscript
# S1 d. MitoCarta share of signal per mouse.

source(here::here("R", "panels.R"))

inputs <- c(
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
de <- readRDS(inputs[["contrasts"]])
share <- tibble(sample_id = names(de$pct_mito), pct_mito = de$pct_mito)
panel <- per_mouse_plot(
  share, as_tibble(de$targets), c(pct_mito = "MitoCarta share"),
  "Pellet mitochondrial share", "About 58% in males, 47% in females; no treatment shift",
  "% of signal"
) +
  theme(strip.text = element_blank())
panel_data <- share
panel_note <- "MitoCarta share of signal per mouse."
save_panel(panel, "F01", "supp", "S1_D_mito_share", 70, 65)
