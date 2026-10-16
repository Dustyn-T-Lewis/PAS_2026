# S7 b: every male candidate shown, which routes it passed, its log2 FC and its STRING
# neighbourhood (06_Candidates).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(candidates = output_path("02_Statistics", "06_Candidates", "candidates.rds"))
cand <- readRDS(inputs[["candidates"]])
routes <- candidate_routes(cand, "male", "Male candidates by route")

save_panel(routes$plot, "S7_Candidates", "b_male_routes", 90, 165)
list(plot = routes$plot, data = routes$data, inputs = inputs)
