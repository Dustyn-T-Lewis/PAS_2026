#!/usr/bin/env Rscript
# S2 d. MDS of the mitochondrial universe (limpa, quantification standard errors).

source(here::here("R", "panels.R"))

inputs <- c(
  proteins = here("01_Preprocess", "02_Quantification", "c_data", "proteins.rds"),
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
de <- readRDS(inputs[["contrasts"]])
proteins <- readRDS(inputs[["proteins"]])$proteins
panel <- mds_plot(
  proteins[de$universes$mito, ], as_tibble(de$targets),
  "Mitochondrial universe ordination"
) +
  labs(subtitle = str_glue("limpa MDS, {sum(de$universes$mito)} proteins")) +
  theme(legend.position = "right")
panel_data <- panel$data
panel_note <- "MDS coordinates, mitochondrial universe (limpa plotMDSUsingSEs)."
save_panel(panel, "F01", "supp", "S2_D_mito_mds", 90, 65)
