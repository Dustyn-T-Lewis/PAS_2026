#!/usr/bin/env Rscript
# How well each gene-set score alone separates the two groups of each comparison: AUC with DeLong's
# interval and the exact Wilcoxon p on the singscores, BH within comparison. Descriptive.

source(here::here("R", "helpers.R"))

out <- stage_paths("03_Pathway_Enrichment", "03_Classify")
inputs <- c(
  scores = here("03_Pathway_Enrichment", "01_Scores", "c_data", "set_scores.rds"),
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
sc <- readRDS(inputs[["scores"]])
de <- readRDS(inputs[["contrasts"]])
m <- rbind(sc$scores$full, sc$scores$mito)

auc <- auc_screen(m, de$targets$group) |>
  rename(set = feature) |>
  mutate(collection = sc$collection[set], .after = set)
plot_index <- write_pdf(
  file.path(out$reports, "03_classify.pdf"),
  auc_pages(rename(auc, feature = set), m, de$targets$group)
)

write_workbook(
  list(summary = auc_summary(auc), auc = auc, plot_index = plot_index),
  c(
    "Per comparison: sets, p ≤ 0.05, FDR ≤ 0.10, perfect separations and the number chance gives.",
    "Every set and comparison: AUC, DeLong 95% interval, exact Wilcoxon p, BH FDR.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)
auc_summary(auc)
sessionInfo()
