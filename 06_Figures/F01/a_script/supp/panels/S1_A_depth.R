#!/usr/bin/env Rscript
# S1 a. Precursors and protein groups quantified per mouse after filtering.

source(here::here("R", "panels.R"))

inputs <- c(
  filtering = here("01_Preprocess", "01_Filtering", "c_data", "01_filter.xlsx"),
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
targets <- as_tibble(readRDS(inputs[["contrasts"]])$targets)
depth <- readxl::read_excel(inputs[["filtering"]], "run_qc")
panel <- per_mouse_plot(
  depth, targets, c(precursors = "Precursors", proteins = "Protein groups"),
  "Data depth", "After filtering; females carry more", "Count"
)
panel_data <- select(depth, sample_id, group, precursors, proteins)
panel_note <- "Precursors and protein groups per mouse (01_Preprocess/01_Filtering)."
save_panel(panel, "F01", "supp", "S1_A_depth", 90, 65)
