# S3 b: male PAS vs vehicle. Above, every protein in the full proteome as log2 FC against -log10 p
# (01_Fit); below, the proteins at Π <= 0.05 across every mouse, centred on each sex's vehicle
# mean, in one slice per co-expression cluster (02_Enrichment) with its direction on the left,
# grey where the protein was not detected in that mouse (02_Quantification), and each cluster's
# ORA bars (02_Enrichment).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))
source(here("04_Figures", "F02_Proteins", "a_script", "supp", "s3_heatmap.R"))

inputs <- c(
  fit = output_path("02_Statistics", "01_Fit", "fit.rds"),
  proteins = output_path("01_Processing", "02_Quantification", "proteins.rds"),
  enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds")
)
panel <- pi_heatmap_panel(
  readRDS(inputs[["fit"]]), readRDS(inputs[["proteins"]]),
  readRDS(inputs[["enrichment"]]), SEX_CONTRASTS[["male"]]
)

save_figure(panel$plot, "F02_Proteins", "S3b_pi_male", panel$width_mm, panel$height_mm)
list(plot = panel$plot, data = panel$data)
