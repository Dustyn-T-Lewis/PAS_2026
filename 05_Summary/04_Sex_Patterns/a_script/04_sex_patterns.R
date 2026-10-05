#!/usr/bin/env Rscript
# Which responses to PAS are shared, which are detected in one sex, and which differ by sex, for
# proteins, pathways and modules. Labels come from each sex's effect and the direct interaction
# test (treatment_by_sex): only a significant interaction says the sexes differ; a hit in one sex
# alone is "detected in" that sex. Nominal p ≤ 0.05; every row carries its FDRs.

source(here::here("R", "helpers.R"))

out <- stage_paths("05_Summary", "04_Sex_Patterns")
xl <- \(path, sheet) readxl::read_excel(here(path), sheet, guess_max = 1e5)
inputs <- c(
  proteins = here("02_Differential_Expression/02_Contrasts/c_data/02_contrasts.xlsx"),
  pathways = here("03_Pathway_Enrichment/02_Contrasts/c_data/02_contrasts.xlsx"),
  modules = here("04_Network/02_Contrasts/c_data/02_contrasts.xlsx")
)

wide <- function(x, id) {
  x |>
    filter(contrast %in% PRIMARY) |>
    mutate(contrast = recode(contrast,
      PAS_vs_VEH_female = "f", PAS_vs_VEH_male = "m",
      treatment_by_sex = "int"
    )) |>
    pivot_wider(id_cols = all_of(id), names_from = contrast, values_from = c(effect, p, fdr)) |>
    mutate(label = sex_label(effect_f, p_f, effect_m, p_m, p_int))
}
proteins <- xl(inputs[["proteins"]], "results") |>
  transmute(
    level = "protein", universe, set = NA_character_, feature = protein, gene, contrast,
    effect = logFC, p = P.Value, fdr = adj.P.Val
  ) |>
  wide(c("level", "universe", "feature", "gene"))
pathways <- xl(inputs[["pathways"]], "fgsea") |>
  transmute(
    level = "pathway", universe, collection, feature = term, contrast, effect = score, p,
    fdr = padj
  ) |>
  wide(c("level", "universe", "collection", "feature"))
modules <- xl(inputs[["modules"]], "results") |>
  transmute(
    level = "module", feature = module, contrast, effect = estimate, p = P.Value,
    fdr = adj.P.Val
  ) |>
  wide(c("level", "feature"))

patterns <- bind_rows(proteins, pathways, modules) |>
  relocate(level, universe, collection, feature, gene, label)
LABELS <- c(
  "sex-differential", "shared", "detected in females", "detected in males",
  "opposite, interaction not significant", "neither"
)
summary <- patterns |>
  count(level, universe, collection, label) |>
  mutate(label = factor(label, LABELS)) |>
  arrange(level, universe, collection, label) |>
  pivot_wider(names_from = label, values_from = n, values_fill = 0)

colours <- set_names(c("#7B3294", "#1B7837", "#B2182B", "#2166AC", "#E08214", "grey80"), LABELS)
scatter <- function(x, title) {
  \() {
    par(mar = c(4, 4, 3, 1))
    x <- arrange(x, label != "neither")
    plot(x$effect_m, x$effect_f,
      pch = 19, cex = 0.5, col = colours[x$label],
      xlab = "effect in males (PAS vs VEH)", ylab = "effect in females (PAS vs VEH)",
      main = title, cex.main = 0.85
    )
    abline(h = 0, v = 0, lty = 3)
    abline(0, 1, col = "grey50")
    legend("topleft", LABELS, pch = 19, col = colours, bty = "n", cex = 0.6)
  }
}
groups <- distinct(patterns, level, universe, collection)
pages <- pmap(groups, \(level, universe, collection) {
  x <- filter(patterns, level == !!level, universe %in% !!universe, collection %in% !!collection)
  scatter(x, paste(na.omit(c(level, universe, collection)), collapse = ", "))
})
names(pages) <- paste(
  "female against male effect,", groups$level, coalesce(groups$universe, ""),
  coalesce(groups$collection, "")
)
plot_index <- write_pdf(file.path(out$reports, "04_sex_patterns.pdf"), pages, width = 8, height = 8)

write_workbook(
  list(summary = summary, patterns = patterns, plot_index = plot_index),
  c(
    "Counts of each label per level, universe and collection.",
    "Every protein, pathway, module: female, male, interaction effect, p, FDR, label.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)
summary
sessionInfo()
