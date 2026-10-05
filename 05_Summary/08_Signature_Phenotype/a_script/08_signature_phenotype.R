#!/usr/bin/env Rscript
# Whether each sex's PAS signature tracks muscle phenotype. Three parts: (1) the PAS effect on each
# phenotype within each sex, in SD units of all twenty mice (Welch t, 95% CI), with the
# treatment-by-sex term from lm; (2) a signature score per mouse, singscore::simpleScore on the
# signature proteins (p ≤ 0.05 for PAS in that sex, up and down sets, 05_Summary/07_Signature);
# (3) each score against each phenotype with the models the protein associations use: pooled
# within-group slope (group means in the design), females only, and the male-minus-female slope
# difference. The signatures were chosen on these mice; keeping group means in every model removes
# the PAS mean difference, but the slopes stay exploratory.

source(here::here("R", "helpers.R"))

SEXES <- c(PAS_vs_VEH_female = "female", PAS_vs_VEH_male = "male")

out <- stage_paths("05_Summary", "08_Signature_Phenotype")
inputs <- c(
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds"),
  signature = here("05_Summary", "07_Signature", "c_data", "signature.rds"),
  phenotypes = here("00_Input", "phenotypes.xlsx")
)
de <- readRDS(inputs[["contrasts"]])
signature <- readRDS(inputs[["signature"]])
samples <- as_tibble(de$targets) |> select(sample_id, treatment, sex, group)
ph <- read_phenotypes(samples$sample_id)
male <- as.numeric(samples$sex == "M")
female <- samples$sex == "F"

phenotype_effects <- map(PHENOTYPES, \(p) {
  z <- as.numeric(scale(ph[[p]]))
  interaction_p <- summary(lm(z ~ treatment * sex, data = samples))$coefficients[
    "treatmentVEH:sexM", "Pr(>|t|)"
  ]
  map(c(F = "F", M = "M"), \(s) {
    test <- t.test(
      z[samples$sex == s & samples$treatment == "PAS"],
      z[samples$sex == s & samples$treatment == "VEH"]
    )
    tibble(
      phenotype = p, sex = s, difference = unname(diff(rev(test$estimate))),
      low = test$conf.int[1], high = test$conf.int[2], p = test$p.value,
      interaction_p = interaction_p
    )
  }) |>
    list_rbind()
}) |>
  list_rbind()

ranks <- singscore::rankGenes(de$E_norm$full)
scores <- imap(SEXES, \(sex, cn) {
  sets <- filter(signature$ranks, source == cn) |> split(~direction)
  singscore::simpleScore(ranks, upSet = sets$up$protein, downSet = sets$down$protein) |>
    as_tibble(rownames = "sample_id") |>
    transmute(sample_id, signature = sex, score = TotalScore)
}) |>
  list_rbind() |>
  left_join(samples, by = "sample_id")

slope <- function(y, design, coef) {
  fit <- summary(lm(y ~ 0 + design))$coefficients
  row <- paste0("design", coef)
  tibble(estimate = fit[row, 1], se = fit[row, 2], p = fit[row, 4])
}
slopes <- expand_grid(signature = unname(SEXES), phenotype = PHENOTYPES) |>
  pmap(\(signature, phenotype) {
    y <- scores$score[scores$signature == signature]
    x <- ph[[phenotype]]
    bind_rows(
      slope(y, with_phenotype(de$design, x), "phenotype") |> mutate(model = "pooled slope"),
      slope(
        y[female], with_phenotype(de$design[female, c("PAS_F", "VEH_F")], x[female]),
        "phenotype"
      ) |>
        mutate(model = "females only"),
      slope(y, with_sex_slope(de$design, x, male), "phenotype_male_minus_female") |>
        mutate(model = "slope difference, male minus female")
    ) |>
      mutate(signature = signature, phenotype = phenotype, .before = 1)
  }) |>
  list_rbind()

saveRDS(
  list(phenotype_effects = phenotype_effects, scores = scores, slopes = slopes),
  file.path(out$data, "signature_phenotype.rds")
)
write_workbook(
  list(phenotype_effects = phenotype_effects, scores = scores, slopes = slopes),
  c(
    "PAS minus vehicle per phenotype and sex, SD units, Welch 95% CI; lm treatment-by-sex p.",
    "Signature score per mouse (singscore, up and down sets) for each sex's PAS signature.",
    "Score against phenotype: pooled within-group slope, females only, slope difference."
  ),
  out$workbook, inputs
)

phenotype_effects
filter(slopes, p <= ALPHA)
sessionInfo()
