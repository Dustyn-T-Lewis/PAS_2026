# F02 b: male candidate proteins (two or more of the six routes in 06_Candidates) as a STRING
# network, shaded by STRING community and named by each community's best over-represented term.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(candidates = output_path("02_Statistics", "06_Candidates", "candidates.rds"))
cand <- readRDS(inputs[["candidates"]])
n_nodes <- sum(cand$candidates$sex == "male" & cand$candidates$routes >= 2)

plot <- candidate_network(cand$candidates, cand$edges, "male",
  title = "Males: candidate proteins",
  subtitle = str_glue("{n_nodes} proteins passing two or more routes; edges: STRING score ≥ 700")
)

save_panel(plot, "F02_Proteins", "b_male_network", 90, 110)
list(
  plot = plot,
  data = filter(cand$candidates, sex == "male", routes >= 2) |>
    select(gene, logFC, P.Value, lfsr, routes, statistics:pathway, community, name),
  inputs = inputs
)
