#!/usr/bin/env Rscript
# How far PAS moves each mouse's whole proteome, by sex. Each mouse's Euclidean distance from its
# own sex's vehicle centroid on the normalised full universe; a vehicle mouse is measured against
# the centroid of the other four, so it is not compared with itself. RRPP (residual randomisation,
# 999 permutations) then asks whether PAS mice sit further out than vehicle mice and whether that
# shift differs by sex (the sex x treatment term); vehicle distances are the noise floor.

source(here::here("R", "helpers.R"))

PERMUTATIONS <- 999

out <- stage_paths("05_Summary", "06_Magnitude")
inputs <- c(
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
de <- readRDS(inputs[["contrasts"]])
samples <- as_tibble(de$targets) |> select(sample_id, treatment, sex, group)
x <- de$E_norm$full[, samples$sample_id]

distance_from_vehicle <- function(id, sex_of, treatment_of) {
  vehicle <- samples$sample_id[samples$sex == sex_of & samples$treatment == "VEH"]
  reference <- setdiff(vehicle, id)
  sqrt(sum((x[, id] - rowMeans(x[, reference, drop = FALSE]))^2))
}
magnitude <- samples |>
  mutate(distance = pmap_dbl(list(sample_id, sex, treatment), distance_from_vehicle))

fit <- RRPP::lm.rrpp(distance ~ sex * treatment,
  data = as.data.frame(magnitude), iter = PERMUTATIONS, seed = SEED, print.progress = FALSE
)
tests <- anova(fit)$table |>
  as.data.frame() |>
  rownames_to_column("term") |>
  as_tibble() |>
  select(term, df = Df, r2 = Rsq, f = `F`, p = `Pr(>F)`)
medians <- summarise(magnitude, median_distance = median(distance), .by = c(sex, treatment))

saveRDS(list(magnitude = magnitude, tests = tests), file.path(out$data, "magnitude.rds"))
write_workbook(
  list(magnitude = magnitude, tests = tests, medians = medians),
  c(
    "Distance of each mouse from its own sex's vehicle centroid (leave-one-out for vehicle mice).",
    "RRPP ANOVA on distance: sex, treatment and their interaction (999 permutations).",
    "Median distance per sex and treatment."
  ),
  out$workbook, inputs
)

tests
medians
sessionInfo()
