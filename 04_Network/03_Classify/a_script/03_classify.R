#!/usr/bin/env Rscript
# How well each module eigengene alone separates the two groups of each comparison: AUC with
# DeLong's interval and the exact Wilcoxon p, BH within comparison. Descriptive.

source(here::here("R", "helpers.R"))

out <- stage_paths("04_Network", "03_Classify")
inputs <- c(
  modules = here("04_Network", "01_Build_Modules", "c_data", "modules.rds"),
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
mods <- readRDS(inputs[["modules"]])
de <- readRDS(inputs[["contrasts"]])
m <- t(as.matrix(mods$eigengenes))

auc <- auc_screen(m, de$targets$group) |> rename(module = feature)
plot_index <- write_pdf(
  file.path(out$reports, "03_classify.pdf"),
  auc_pages(rename(auc, feature = module), m, de$targets$group, n = 9)
)

write_workbook(
  list(summary = auc_summary(auc), auc = auc, plot_index = plot_index),
  c(
    "Per comparison: p ≤ 0.05, FDR ≤ 0.10, perfect separations, chance's count.",
    "Every module and comparison: AUC, DeLong 95% interval, exact Wilcoxon p, BH FDR.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)
auc_summary(auc)
sessionInfo()
