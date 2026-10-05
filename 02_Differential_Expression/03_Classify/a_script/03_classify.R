#!/usr/bin/env Rscript
# How well each protein alone separates the two groups of each comparison: AUC with DeLong's
# interval and the exact Wilcoxon p, BH within comparison. Descriptive: at five against five many
# proteins separate perfectly by chance, and the summary says how many to expect.

source(here::here("R", "helpers.R"))

out <- stage_paths("02_Differential_Expression", "03_Classify")
inputs <- c(
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
de <- readRDS(inputs[["contrasts"]])
m <- de$E_norm$full
gene <- set_names(str_split_i(de$genes$Genes, ";", 1), rownames(de$genes))

auc <- auc_screen(m, de$targets$group) |>
  rename(protein = feature) |>
  mutate(gene = gene[protein], .after = protein)
plot_index <- write_pdf(
  file.path(out$reports, "03_classify.pdf"),
  auc_pages(rename(auc, feature = protein), m, de$targets$group, \(f) gene[f])
)

write_workbook(
  list(summary = auc_summary(auc), auc = auc, plot_index = plot_index),
  c(
    "Per comparison: p ≤ 0.05, FDR ≤ 0.10, perfect separations, chance's count.",
    "Every protein and comparison: AUC, DeLong 95% interval, exact Wilcoxon p, BH FDR.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)
auc_summary(auc)
sessionInfo()
