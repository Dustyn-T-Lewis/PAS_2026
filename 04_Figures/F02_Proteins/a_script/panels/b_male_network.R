# F02 b: male candidate proteins (06_Candidates) as a STRING network, shaded by STRING community
# and labelled by each community's ORA term.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(candidates = output_path("02_Statistics", "06_Candidates", "candidates.rds"))
cand <- readRDS(inputs[["candidates"]])
network <- candidate_network(cand, "male", "Males: candidate proteins")

save_panel(network$plot, "F02_Proteins", "b_male_network", 90, 110)
list(plot = network$plot, data = network$data, inputs = inputs)
