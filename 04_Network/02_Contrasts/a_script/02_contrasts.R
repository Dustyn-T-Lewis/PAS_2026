#!/usr/bin/env Rscript
# Whether each module's eigengene moved with treatment or sex: an ordinary linear model on the four
# group means per module, the five contrasts by emmeans, BH across the nine modules within each
# contrast (grey has no eigengene). Nine modules are too few for limma's empirical Bayes prior.
# WGCNA's module-trait correlations for treatment and sex sit beside them.

source(here::here("R", "helpers.R"))
library(emmeans)
library(WGCNA)

out <- stage_paths("04_Network", "02_Contrasts")
inputs <- c(
  modules = here("04_Network", "01_Build_Modules", "c_data", "modules.rds"),
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
mods <- readRDS(inputs[["modules"]])
de <- readRDS(inputs[["contrasts"]])
eigengenes <- mods$eigengenes
group <- factor(de$targets$group, levels = rownames(de$contrasts))
stopifnot(identical(rownames(eigengenes), rownames(de$targets)))

contrast_list <- map(CONTRASTS, \(cn) unname(de$contrasts[, cn])) |> set_names(CONTRASTS)
results <- map(colnames(eigengenes), \(m) {
  fit <- lm(eigengenes[[m]] ~ 0 + group)
  contrast(emmeans(fit, ~group), method = contrast_list) |>
    as_tibble() |>
    transmute(
      module = m, contrast = as.character(contrast), estimate, SE, df, t = t.ratio,
      P.Value = p.value
    )
}) |>
  list_rbind() |>
  mutate(adj.P.Val = p.adjust(P.Value, "BH"), .by = contrast) |>
  left_join(select(mods$sizes, module, proteins, mito_share), by = "module")

summary <- results |>
  reframe(tier_counts(P.Value, adj.P.Val), .by = contrast)

traits <- cbind(
  treatment_PAS = as.numeric(de$targets$treatment == "PAS"),
  sex_male = as.numeric(de$targets$sex == "M")
)
module_trait <- expand_grid(module = colnames(eigengenes), trait = colnames(traits)) |>
  mutate(
    r = map2_dbl(module, trait, \(m, tr) cor(eigengenes[[m]], traits[, tr])),
    P.Value = corPvalueStudent(r, nrow(traits))
  )

r_matrix <- pivot_wider(module_trait, id_cols = module, names_from = trait, values_from = r) |>
  tibble::column_to_rownames("module") |>
  as.matrix()
p_matrix <- module_trait |>
  pivot_wider(id_cols = module, names_from = trait, values_from = P.Value) |>
  tibble::column_to_rownames("module") |>
  as.matrix()
overview <- write_pdf(file.path(out$reports, "02_contrasts_module_trait.pdf"), list(
  "module-trait correlation, treatment and sex" = \() {
    labels <- matrix(sprintf("%.2f\n(p %.3f)", r_matrix, p_matrix), nrow(r_matrix))
    par(mar = c(6, 9, 3, 3))
    labeledHeatmap(r_matrix,
      xLabels = colnames(r_matrix), yLabels = rownames(r_matrix),
      colors = blueWhiteRed(50), textMatrix = labels, zlim = c(-1, 1), cex.text = 0.7,
      main = "Module-trait correlations (Pearson r, corPvalueStudent)"
    )
  }
)) |>
  mutate(contrast = "all", .before = 1)

per_contrast <- map(CONTRASTS, \(cn) {
  x <- filter(results, contrast == cn)
  write_pdf(file.path(out$reports, str_glue("02_contrasts_{cn}.pdf")), list(
    "eigengenes by group, with this contrast's p and FDR" = \() {
      par(mfrow = c(3, 3), mar = c(3, 4, 3, 1))
      walk(colnames(eigengenes), \(m) {
        r <- x[x$module == m, ]
        stripchart(eigengenes[[m]] ~ factor(de$targets$group, GROUPS),
          vertical = TRUE,
          method = "jitter", pch = 19, cex = 0.8, col = "grey30", ylab = "eigengene",
          main = sprintf("%s: p %.3f, FDR %.2f", m, r$P.Value, r$adj.P.Val), cex.main = 0.85
        )
      })
    }
  )) |>
    mutate(contrast = cn, .before = 1)
}) |>
  list_rbind()

saveRDS(
  list(results = results, module_trait = module_trait),
  file.path(out$data, "module_contrasts.rds")
)
write_workbook(
  list(
    summary = summary, results = results, module_trait = module_trait,
    plot_index = bind_rows(overview, per_contrast)
  ),
  c(
    "Per contrast: modules tested, p ≤ 0.05, chance, FDR ≤ 0.10 and ≤ 0.05.",
    "Every module and contrast: estimate, SE, df, t, p (emmeans), BH FDR across modules.",
    "WGCNA module-trait r with treatment (PAS = 1) and sex (male = 1), corPvalueStudent.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

filter(results, P.Value <= ALPHA)

sessionInfo()
