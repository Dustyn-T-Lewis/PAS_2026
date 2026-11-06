# What shapes the whole proteome and the phenotypes. Proteome: PCA of the twenty mice on the
# normalised full universe (centred, not scaled), PERMANOVA (vegan::adonis2, Euclidean, sequential
# terms) for sex, treatment and their interaction, a dispersion check (betadisper by group), the R²
# of each design factor, the pellet's MitoCarta share and run depth on PC1 to PC4, and the
# MitoCarta share fitted onto PC1-PC2 (vegan::envfit). Phenotypes: PCA of the six traits, scaled,
# with the same PERMANOVA.

source(here::here("00_Setup", "Shared_functions", "helpers.R"))

N_COMPONENTS <- 4
PERMUTATIONS <- 999

out <- stage_paths("02_Statistics", "04_PCA")
inputs <- c(
  fit = output_path("02_Statistics", "01_Fit", "fit.rds"),
  filtering = output_path("01_Processing", "01_Filtering", "precursors_filtered.rds")
)
fit <- readRDS(inputs[["fit"]])
run_qc <- readRDS(inputs[["filtering"]])$run_qc
samples <- as_tibble(fit$targets) |>
  select(sample_id, treatment, sex, group) |>
  mutate(pct_mito = fit$pct_mito[sample_id], run_depth = run_depth_log2(run_qc, sample_id)) |>
  left_join(
    select(run_qc, sample_id, precursors, protein_groups = proteins, raw_signal_rel = raw_rel),
    by = "sample_id"
  )

# Technical variables per group, with the two-way model each is read against. They are reported,
# not adjusted for: run order was randomised and within-sex fold changes cancel what PAS and
# vehicle share.
TECHNICAL <- c("precursors", "protein_groups", "raw_signal_rel", "pct_mito")
technical_groups <- samples |>
  summarise(across(all_of(TECHNICAL), median), .by = group)
technical_tests <- map(TECHNICAL, \(v) {
  coefs <- summary(lm(samples[[v]] ~ treatment * sex, samples))$coefficients
  tibble(variable = v, term = rownames(coefs)[-1], p = coefs[-1, "Pr(>|t|)"])
}) |>
  list_rbind()

ordinate <- function(x, scale) {
  pca <- prcomp(x, center = TRUE, scale. = scale)
  distances <- dist(if (scale) scale(x) else x)
  set.seed(SEED)
  permanova <- vegan::adonis2(distances ~ sex * treatment,
    data = samples, permutations = PERMUTATIONS, by = "terms"
  ) |>
    as.data.frame() |>
    rownames_to_column("term") |>
    as_tibble() |>
    select(term, df = Df, r2 = R2, f = `F`, p = `Pr(>F)`)
  dispersion <- vegan::betadisper(distances, samples$group)
  set.seed(SEED)
  dispersion_test <- vegan::permutest(dispersion, permutations = PERMUTATIONS)
  list(
    pca = pca, variance = pca$sdev^2 / sum(pca$sdev^2), permanova = permanova,
    dispersion = dispersion, dispersion_p = dispersion_test$tab[1, "Pr(>F)"]
  )
}
proteome <- ordinate(t(fit$E_norm$full[, samples$sample_id]), scale = FALSE)
phenotype <- ordinate(as.matrix(fit$phenotypes[PHENOTYPES]), scale = TRUE)

set.seed(SEED)
mito_fit <- vegan::envfit(proteome$pca$x[, 1:2], samples["pct_mito"], permutations = PERMUTATIONS)
mito_vector <- tibble(
  PC1 = mito_fit$vectors$arrows[1, 1], PC2 = mito_fit$vectors$arrows[1, 2],
  r2 = unname(mito_fit$vectors$r), p = unname(mito_fit$vectors$pvals)
)

codes <- transmute(samples,
  Sex = if_else(sex == "M", 1, -1),
  Treatment = if_else(treatment == "PAS", 1, -1),
  `Treatment × sex` = Treatment * Sex,
  `MitoCarta share` = pct_mito,
  `Run depth (log2 precursors)` = run_depth
)
component_r2 <- expand_grid(component = seq_len(N_COMPONENTS), factor = names(codes)) |>
  mutate(
    r2 = map2_dbl(component, factor, \(k, f) cor(proteome$pca$x[, k], codes[[f]])^2),
    pc = sprintf("PC%d (%.1f%%)", component, 100 * proteome$variance[component])
  )
scores <- function(o, n) {
  as_tibble(o$pca$x[, seq_len(n)]) |>
    bind_cols(samples, .name_repair = "minimal")
}
proteome_scores <- scores(proteome, N_COMPONENTS)
phenotype_scores <- select(scores(phenotype, 2), -c(run_depth, protein_groups, raw_signal_rel))
phenotype_loadings <- as_tibble(phenotype$pca$rotation[, 1:2], rownames = "phenotype")

group <- factor(samples$group, FIGURE_GROUPS)
plot_index <- preview_pdf(file.path(out$reports, "04_pca.pdf"), list(
  "proteome PCA" = \() {
    xy <- proteome$pca$x[, 1:2]
    vegan::ordiplot(xy, type = "n", main = "Proteome PCA")
    points(xy, pch = SEX_SHAPES[samples$sex], col = GROUP_COLOURS[as.character(group)])
    vegan::ordiellipse(xy, group, kind = "sd", col = GROUP_COLOURS[levels(group)])
    plot(mito_fit, col = "grey20")
    legend("topleft", levels(group), col = GROUP_COLOURS[levels(group)], pch = 16, bty = "n")
  },
  "phenotype PCA" = \() biplot(phenotype$pca, xlabs = samples$group, cex = 0.6),
  "proteome dispersion by group" = \() plot(proteome$dispersion, main = "Proteome dispersion")
))

saveRDS(
  list(
    proteome_scores = proteome_scores, proteome_variance = proteome$variance,
    permanova = proteome$permanova, component_r2 = component_r2, mito_vector = mito_vector,
    phenotype_scores = phenotype_scores, phenotype_variance = phenotype$variance,
    phenotype_permanova = phenotype$permanova, phenotype_loadings = phenotype_loadings
  ),
  file.path(out$data, "pca.rds")
)
write_workbook(
  list(
    proteome_scores = mutate(proteome_scores, sex = SEX_NAMES[sex]),
    proteome_variance = tibble(component = seq_along(proteome$variance), share = proteome$variance),
    permanova = mutate(proteome$permanova, dispersion_p = c(
      proteome$dispersion_p, rep(NA, nrow(proteome$permanova) - 1)
    )),
    component_r2 = select(component_r2, component, factor, r2),
    mito_vector = mito_vector,
    technical_groups = technical_groups, technical_tests = technical_tests,
    phenotype_scores = mutate(phenotype_scores, sex = SEX_NAMES[sex]),
    phenotype_loadings = phenotype_loadings,
    phenotype_permanova = mutate(phenotype$permanova, dispersion_p = c(
      phenotype$dispersion_p, rep(NA, nrow(phenotype$permanova) - 1)
    )),
    plot_index = plot_index
  ),
  c(
    "Proteome PCA scores of the twenty mice, normalised full universe, centred, not scaled.",
    "Share of proteome variance per component.",
    str_glue(
      "Proteome PERMANOVA (adonis2, Euclidean, {PERMUTATIONS} permutations, sequential); ",
      "betadisper p, row 1."
    ),
    "R² of each design factor, the MitoCarta share and run depth on PC1 to PC4.",
    "MitoCarta share fitted onto PC1-PC2 (envfit): unit direction, R² and permutation p.",
    "Median run depth, protein groups, relative raw signal and MitoCarta share per group.",
    "Each technical variable against treatment, sex and their interaction (lm p).",
    "Phenotype PCA scores, six traits scaled.",
    "Phenotype PCA loadings on PC1 and PC2.",
    "Phenotype PERMANOVA on scaled traits; betadisper p on row 1.",
    "Which PDF page shows what."
  ),
  out$workbook
)

proteome$permanova
phenotype$permanova
