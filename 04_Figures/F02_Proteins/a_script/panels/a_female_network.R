# F02 a: female candidate proteins (06_Candidates) as a STRING network, shaded by STRING community
# and labelled by each community's ORA term.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(candidates = output_path("02_Statistics", "06_Candidates", "candidates.rds"))
cand <- readRDS(inputs[["candidates"]])
network <- candidate_network(cand, "female", "Females: candidate proteins")

save_panel(network$plot, "F02_Proteins", "a_female_network", 90, 90)
list(plot = network$plot, data = network$data, inputs = inputs)
