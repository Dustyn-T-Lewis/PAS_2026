#!/usr/bin/env Rscript
# Whether one sex's PAS signature moves the same way in the other sex. The proteins at p ≤ 0.05 for
# PAS in one sex, split by direction, are tested as sets in the other sex's PAS contrast with
# limma::fry (self-contained, precision weights, the stored design and contrasts) on the normalised
# full universe; ranks in the other sex's t ordering are kept for barcode plots. A signature that
# moves the same way points to a shared response; one that does not, to a sex-specific one, which
# only the treatment-by-sex contrast tests directly.

source(here::here("R", "helpers.R"))

SEXES <- c(PAS_vs_VEH_female = "female", PAS_vs_VEH_male = "male")

out <- stage_paths("05_Summary", "07_Signature")
inputs <- c(
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
de <- readRDS(inputs[["contrasts"]])
full <- filter(de$results, universe == "full")

signatures <- full |>
  filter(contrast %in% names(SEXES), P.Value <= ALPHA) |>
  mutate(direction = if_else(logFC > 0, "up", "down")) |>
  select(source = contrast, direction, protein)
pairs <- tibble(source = names(SEXES), target = rev(names(SEXES)))

tests <- pmap(pairs, \(source, target) {
  sets <- signatures |>
    filter(source == !!source) |>
    split(~direction) |>
    map(\(x) match(x$protein, rownames(de$E_norm$full)))
  limma::fry(de$E_norm$full, sets,
    design = de$design, contrast = de$contrasts[, target], weights = de$weights$full
  ) |>
    as_tibble(rownames = "direction") |>
    transmute(
      source, target, direction,
      proteins = NGenes, moves = Direction, p = PValue
    )
}) |>
  list_rbind()

ranks <- pmap(pairs, \(source, target) {
  ranking <- filter(full, contrast == target) |>
    transmute(protein, rank = rank(-t) / n())
  signatures |>
    filter(source == !!source) |>
    left_join(ranking, by = "protein") |>
    mutate(target = target, .after = source)
}) |>
  list_rbind()
agreement <- ranks |>
  summarise(
    same_way = mean(if_else(direction == "up", rank < 0.5, rank > 0.5)),
    .by = c(source, target, direction)
  )
tests <- left_join(tests, agreement, by = c("source", "target", "direction"))

saveRDS(list(tests = tests, ranks = ranks), file.path(out$data, "signature.rds"))
write_workbook(
  list(tests = tests, ranks = ranks),
  c(
    "fry test of each sex's PAS signature (p ≤ 0.05, by direction) in the other sex's contrast.",
    "Each signature protein's rank (0 = most up, 1 = most down) in the other sex's t ordering."
  ),
  out$workbook, inputs
)

tests
sessionInfo()
