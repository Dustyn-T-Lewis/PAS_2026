# F02 d: every mouse for the male candidates (06_Candidates), as log2 abundance against the
# own-sex vehicle mean in all four groups; open circles mark mice in which the protein was not
# detected.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))
source(here("04_Figures", "F02_Proteins", "a_script", "candidates.R"))

inputs <- c(
  candidates = output_path("02_Statistics", "06_Candidates", "candidates.rds"),
  fit = output_path("02_Statistics", "01_Fit", "fit.rds"),
  proteins = output_path("01_Processing", "02_Quantification", "proteins.rds")
)
dots <- candidate_dots(
  readRDS(inputs[["candidates"]])$candidates, readRDS(inputs[["fit"]]),
  readRDS(inputs[["proteins"]]), "male", "Male candidates, every mouse"
)

save_figure(dots$plot, "F02_Proteins", "d_male_candidates", 180, 62)
list(plot = dots$plot, data = dots$data)
