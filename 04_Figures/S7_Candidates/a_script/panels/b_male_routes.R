# S7 b: every male candidate passing two or more routes, which routes it passed, its log2 FC and
# its STRING neighbourhood.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(candidates = output_path("02_Statistics", "06_Candidates", "candidates.rds"))
cand <- readRDS(inputs[["candidates"]])$candidates

plot <- candidate_routes(cand, "male",
  title = "Male candidates by route",
  subtitle = "Filled: route passed; right: STRING neighbourhood; sorted by routes, then p"
)

save_panel(plot, "S7_Candidates", "b_male_routes", 90, 165)
list(
  plot = plot,
  data = cand |>
    filter(sex == "male", routes >= 2) |>
    arrange(desc(routes), P.Value) |>
    select(gene, routes, all_of(names(CANDIDATE_ROUTES)), logFC, P.Value, lfsr, name),
  inputs = inputs
)
