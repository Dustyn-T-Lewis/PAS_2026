#!/usr/bin/env Rscript
# S1 b. Plasma and contaminant share of signal per mouse, before removal.

source(here::here("R", "panels.R"))

inputs <- c(
  filtering = here("01_Preprocess", "01_Filtering", "c_data", "01_filter.xlsx"),
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
targets <- as_tibble(readRDS(inputs[["contrasts"]])$targets)
contamination <- readxl::read_excel(inputs[["filtering"]], "contaminant_per_run")
panel <- per_mouse_plot(
  contamination, targets,
  c(pct_plasma = "Plasma", pct_contaminant = "Contaminants"),
  "Contamination", "Share of signal before removal; no pattern by group", "% of signal"
)
panel_data <- contamination
panel_note <- "Plasma and contaminant share of signal per mouse before removal (01_Filtering)."
save_panel(panel, "F01", "supp", "S1_B_contamination", 90, 65)
