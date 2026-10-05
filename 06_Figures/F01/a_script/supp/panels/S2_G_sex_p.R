#!/usr/bin/env Rscript
# S2 g. p-value histograms: baseline sex p-values.

source(here::here("R", "panels.R"))

LABELS <- c(
  sex_in_VEH = "Sex in VEH\n(male − female)", sex_in_PAS = "Sex in PAS\n(male − female)"
)
inputs <- c(
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds"),
  signal = here("05_Summary", "01_Signal", "c_data", "01_signal.xlsx")
)
signal <- readxl::read_excel(inputs[["signal"]], "signal") |>
  filter(level == "protein", question == "contrast", universe == "full")
results <- filter(readRDS(inputs[["contrasts"]])$results, universe == "full")
panel <- p_histogram(results, LABELS, signal, "Baseline sex p-values") +
  labs(subtitle = "Full proteome")
panel_data <- filter(signal, term %in% names(LABELS))
panel_note <- "p ≤ 0.05, chance, FDR and true-effect share per contrast (05_Summary/01_Signal)."
save_panel(panel, "F01", "supp", "S2_G_sex_p", 80, 65)
