#!/usr/bin/env Rscript
# Whether each module eigengene tracks each phenotype, by ordinary least squares under the three
# models (pooled slope, females only, male-minus-female slope difference), BH across the nine
# modules. The pooled module-phenotype correlations of WGCNA's tutorial are kept on a labelled
# descriptive sheet: pooling across groups restates group differences. DGCA is supplementary.

source(here::here("R", "helpers.R"))
library(WGCNA)

out <- stage_paths("04_Network", "04_Associate")
inputs <- c(
  modules = here("04_Network", "01_Build_Modules", "c_data", "modules.rds"),
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds"),
  phenotypes = here("00_Input", "phenotypes.xlsx")
)
mods <- readRDS(inputs[["modules"]])
de <- readRDS(inputs[["contrasts"]])
eigengenes <- mods$eigengenes
ph <- read_phenotypes(rownames(de$targets))
group <- factor(de$targets$group)
male <- as.numeric(de$targets$sex == "M")
female <- de$targets$sex == "F"

slope <- function(fit, term) {
  s <- coef(summary(fit))[term, ]
  tibble(slope_per_sd = s[[1]], se = s[[2]], t = s[[3]], P.Value = s[[4]])
}
association <- expand_grid(module = colnames(eigengenes), phenotype = PHENOTYPES) |>
  pmap(\(module, phenotype) {
    y <- eigengenes[[module]]
    z <- as.numeric(scale(ph[[phenotype]]))
    zf <- as.numeric(scale(ph[[phenotype]][female]))
    sex_slope <- z * male
    bind_rows(
      mutate(slope(lm(y ~ 0 + group + z), "z"), model = "pooled slope"),
      mutate(slope(lm(y[female] ~ 0 + droplevels(group[female]) + zf), "zf"),
        model = "females only"
      ),
      mutate(slope(lm(y ~ 0 + group + z + sex_slope), "sex_slope"),
        model = "slope difference, male minus female"
      )
    ) |>
      mutate(module = module, phenotype = phenotype, .before = 1)
  }) |>
  list_rbind() |>
  mutate(adj.P.Val = p.adjust(P.Value, "BH"), .by = c(phenotype, model))

pooled_correlation <- expand_grid(module = colnames(eigengenes), phenotype = PHENOTYPES) |>
  mutate(
    r = map2_dbl(module, phenotype, \(m, p) cor(eigengenes[[m]], ph[[p]])),
    P.Value = corPvalueStudent(r, nrow(eigengenes)),
    note = "descriptive: pooled across groups, restates group differences"
  )

set.seed(SEED)
dgca <- differential_correlation(t(as.matrix(eigengenes)), ph, adjust = "BH") |>
  rename(module = feature)

summary <- association |>
  reframe(tier_counts(P.Value, adj.P.Val), .by = c(phenotype, model))

t_page <- function(m) {
  x <- filter(association, model == m)
  tm <- pivot_wider(x, id_cols = module, names_from = phenotype, values_from = t) |>
    tibble::column_to_rownames("module") |>
    as.matrix()
  pm <- pivot_wider(x, id_cols = module, names_from = phenotype, values_from = P.Value) |>
    tibble::column_to_rownames("module") |>
    as.matrix()
  \() {
    par(mar = c(9, 8, 3, 3))
    labeledHeatmap(tm / max(abs(tm)),
      xLabels = colnames(tm), yLabels = rownames(tm),
      colors = blueWhiteRed(50), textMatrix = matrix(sprintf("%.1f\n(p %.3f)", tm, pm), nrow(tm)),
      zlim = c(-1, 1), cex.text = 0.6, main = str_glue("Eigengene against phenotype, t ({m})")
    )
  }
}
plot_index <- write_pdf(
  file.path(out$reports, "04_associate.pdf"),
  set_names(
    map(unique(association$model), t_page),
    str_glue("module-phenotype t, {unique(association$model)}")
  )
)

write_workbook(
  list(
    summary = summary, association = association, pooled_correlation = pooled_correlation,
    dgca = dgca, plot_index = plot_index
  ),
  c(
    "Per phenotype and model: modules tested, p ≤ 0.05, chance, FDR ≤ 0.10 and ≤ 0.05.",
    "Every module, phenotype and model: slope per SD, SE, t, p (lm), BH across modules.",
    "Descriptive only: pooled module-phenotype Pearson r with corPvalueStudent.",
    "Supplementary: DGCA per module, comparison and phenotype, BH.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

filter(association, P.Value <= ALPHA)

sessionInfo()
