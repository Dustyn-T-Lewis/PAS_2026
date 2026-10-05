#!/usr/bin/env Rscript
# S2 f. p-value histograms: mitochondrial universe p-values.

source(here::here("R", "panels.R"))

LABELS <- c(
  PAS_vs_VEH_female = "PAS vs VEH\nfemale", PAS_vs_VEH_male = "PAS vs VEH\nmale",
  treatment_by_sex = "Treatment ×\nsex"
)
inputs <- c(
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds"),
  signal = here("05_Summary", "01_Signal", "c_data", "01_signal.xlsx")
)
signal <- readxl::read_excel(inputs[["signal"]], "signal") |>
  filter(level == "protein", question == "contrast", universe == "mito")
results <- filter(readRDS(inputs[["contrasts"]])$results, universe == "mito")
panel <- p_histogram(results, LABELS, signal, "Mitochondrial universe p-values") +
  labs(subtitle = "Dashed: no effect; share of true effects from limma propTrueNull")
panel_data <- filter(signal, term %in% names(LABELS))
panel_note <- "p ≤ 0.05, chance, FDR and true-effect share per contrast (05_Summary/01_Signal)."
save_panel(panel, "F01", "supp", "S2_F_mito_p", 110, 65)
