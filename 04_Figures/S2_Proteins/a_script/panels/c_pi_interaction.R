# S2 c: treatment × sex, the proteins at Π <= 0.05 in the full proteome across every mouse,
# centred on each sex's vehicle mean and split by direction (01_Fit); grey where the protein was
# not detected in that mouse (02_Quantification), with each block's ORA bars (02_Enrichment).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(
  fit = output_path("02_Statistics", "01_Fit", "fit.rds"),
  proteins = output_path("01_Processing", "02_Quantification", "proteins.rds"),
  enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds")
)
panel <- pi_heatmap_panel(
  readRDS(inputs[["fit"]]), readRDS(inputs[["proteins"]]),
  readRDS(inputs[["enrichment"]])$block_ora, "treatment_by_sex"
)

save_panel(panel$plot, "S2_Proteins", "c_pi_interaction", panel$width_mm, panel$height_mm)
list(plot = panel$plot, data = panel$data, inputs = inputs)
