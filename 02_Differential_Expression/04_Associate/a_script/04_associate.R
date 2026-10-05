#!/usr/bin/env Rscript
# Which proteins track each phenotype across mice of the same group, fitted three ways with limpa's
# dpcDE: one slope shared by all 20 mice beside the four group means; the ten females alone; and the
# male-minus-female difference in slope, the direct test of whether the relationship differs by
# sex. DGCA's correlation classes are a supplementary sheet.

source(here::here("R", "helpers.R"))
library(limma)
library(limpa)

out <- stage_paths("02_Differential_Expression", "04_Associate")
inputs <- c(
  proteins = here("01_Preprocess", "02_Quantification", "c_data", "proteins.rds"),
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds"),
  phenotypes = here("00_Input", "phenotypes.xlsx")
)

proteins <- readRDS(inputs[["proteins"]])$proteins
de <- readRDS(inputs[["contrasts"]])
proteins$E <- de$E_norm$full
ph <- read_phenotypes(rownames(de$targets))
male <- as.numeric(de$targets$sex == "M")
female <- de$targets$sex == "F"
female_design <- de$design[female, c("PAS_F", "VEH_F")]

annotation <- as_tibble(de$genes, rownames = "protein") |>
  transmute(protein, gene = str_split_i(Genes, ";", 1), is_mito)

fit_coef <- function(elist, design, coef) {
  eBayes(dpcDE(elist, design, sample.weights = TRUE)) |>
    topTable(coef = coef, number = Inf, sort.by = "none") |>
    rownames_to_column("protein") |>
    transmute(protein, log2fc_per_sd = logFC, t, P.Value, adj.P.Val)
}
association <- map(PHENOTYPES, \(p) {
  bind_rows(
    fit_coef(proteins, with_phenotype(de$design, ph[[p]]), "phenotype") |>
      mutate(model = "pooled slope"),
    fit_coef(proteins[, female], with_phenotype(female_design, ph[[p]][female]), "phenotype") |>
      mutate(model = "females only"),
    fit_coef(proteins, with_sex_slope(de$design, ph[[p]], male), "phenotype_male_minus_female") |>
      mutate(model = "slope difference, male minus female")
  ) |>
    mutate(phenotype = p, .before = 1)
}) |>
  list_rbind() |>
  left_join(annotation, by = "protein") |>
  relocate(phenotype, model, protein, gene, is_mito)

summary <- association |>
  reframe(tier_counts(P.Value, adj.P.Val),
    share_true_effects = 1 - propTrueNull(P.Value),
    .by = c(phenotype, model)
  )
csa <- summary$phenotype == "csa_mean" & summary$model == "pooled slope"
stopifnot(summary$fdr_le_0.10[csa] == 327)

set.seed(SEED)
dgca <- differential_correlation(de$E_norm$full, ph) |>
  rename(protein = feature)

plot_index <- map(PHENOTYPES, \(p) {
  x <- filter(association, phenotype == p)
  d <- filter(dgca, phenotype == p)
  write_pdf(file.path(out$reports, str_glue("04_associate_{p}.pdf")), list(
    "volcanoes: pooled slope, females only, slope difference" = \() {
      par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
      walk(unique(x$model), \(m) {
        y <- filter(x, model == m)
        plot(y$log2fc_per_sd, -log10(y$P.Value),
          pch = 19, cex = 0.4,
          col = if_else(y$adj.P.Val <= 0.10, "#B2182B", "grey70"),
          xlab = "log2 change per SD of phenotype", ylab = "-log10 p", main = m, cex.main = 0.9
        )
        abline(v = 0, h = -log10(ALPHA), lty = 3)
      })
      mtext(str_glue("{p}: red FDR ≤ 0.10"), side = 3, outer = TRUE, line = -1.2, cex = 0.8)
    },
    "DGCA: r in each group of each comparison" = \() {
      par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
      iwalk(COMPARISONS, \(groups, cmp) {
        y <- filter(d, comparison == cmp)
        plot(y$r_second, y$r_first,
          pch = if_else(y$q_diff <= 0.10, 19, 1), cex = 0.4,
          xlim = c(-1, 1), ylim = c(-1, 1), xlab = str_glue("r in {groups[2]}"),
          ylab = str_glue("r in {groups[1]}"), main = cmp, cex.main = 0.9
        )
        abline(h = 0, v = 0, lty = 3)
        abline(0, 1, col = "grey60")
      })
    }
  )) |>
    mutate(phenotype = p, .before = 1)
}) |>
  list_rbind()

saveRDS(list(association = association), file.path(out$data, "protein_association.rds"))
write_workbook(
  list(summary = summary, association = association, dgca = dgca, plot_index = plot_index),
  c(
    "Per phenotype and model: tests, p ≤ 0.05, chance, FDR ≤ 0.10 and ≤ 0.05, true-effect share.",
    "Every protein, phenotype and model: log2 change per SD of phenotype, moderated t, p, BH FDR.",
    "Supplementary: DGCA per protein, comparison, phenotype (r by group, z, p, q, class).",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

summary |> filter(model != "pooled slope" | fdr_le_0.10 > 0)

sessionInfo()
