# S3 c: treatment × sex, the proteins at Π <= 0.05 in the full proteome across every mouse,
# centred on each sex's vehicle mean, split by direction (01_Fit) and co-expression cluster
# (02_Enrichment); grey where the protein was not detected in that mouse (02_Quantification),
# with each cluster's ORA bars (02_Enrichment).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(
  fit = output_path("02_Statistics", "01_Fit", "fit.rds"),
  proteins = output_path("01_Processing", "02_Quantification", "proteins.rds"),
  enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds")
)
panel <- pi_heatmap_panel(
  readRDS(inputs[["fit"]]), readRDS(inputs[["proteins"]]),
  readRDS(inputs[["enrichment"]]), "treatment_by_sex"
)

save_figure(panel$plot, "F02_Proteins", "S3c_pi_interaction", panel$width_mm, panel$height_mm)
list(plot = panel$plot, data = panel$data)
