#!/usr/bin/env Rscript
# What shapes the whole proteome and the phenotypes. Proteome: PCA of the twenty mice on the
# normalised full universe (centred, not scaled), PERMANOVA (vegan::adonis2, Euclidean, 999
# permutations, sequential terms) for sex, treatment and their interaction, a dispersion check
# (betadisper by group), the R² of each design factor, the pellet's MitoCarta share and run depth
# on PC1 to PC4, and the MitoCarta share fitted onto PC1-PC2 (vegan::envfit). Phenotypes: PCA of
# the six traits, scaled, with the same PERMANOVA.

source(here::here("R", "helpers.R"))

N_COMPONENTS <- 4
PERMUTATIONS <- 999

out <- stage_paths("00_Trunk", "10_Global")
inputs <- c(
  fit = here("00_Trunk", "07_Fit", "c_data", "fit.rds"),
  filtering = here("00_Trunk", "01_Filtering", "c_data", "precursors_filtered.rds")
)
fit <- readRDS(inputs[["fit"]])
run_qc <- readRDS(inputs[["filtering"]])$run_qc
samples <- as_tibble(fit$targets) |>
  select(sample_id, treatment, sex, group) |>
  mutate(
    pct_mito = fit$pellet_mito[sample_id],
    precursors = run_qc$precursors[match(sample_id, run_qc$sample_id)],
    protein_groups = run_qc$proteins[match(sample_id, run_qc$sample_id)],
    raw_signal_rel = run_qc$raw_rel[match(sample_id, run_qc$sample_id)]
  )
stopifnot(!anyNA(samples$precursors))

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
  set.seed(SEED)
  dispersion <- vegan::betadisper(distances, samples$group) |>
    vegan::permutest(permutations = PERMUTATIONS)
  list(
    pca = pca, variance = pca$sdev^2 / sum(pca$sdev^2), permanova = permanova,
    dispersion_p = dispersion$tab[1, "Pr(>F)"]
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
  `Sex × treatment` = Sex * Treatment,
  `MitoCarta share` = pct_mito,
  `Run depth (precursors)` = precursors
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
phenotype_scores <- scores(phenotype, 2)
phenotype_loadings <- as_tibble(phenotype$pca$rotation[, 1:2], rownames = "phenotype")

pca_page <- function(s, o, title) {
  ggplot(s, aes(PC1, PC2, colour = factor(group, FIGURE_GROUPS), shape = sex)) +
    geom_point(size = 2) +
    scale_colour_manual(values = GROUP_COLOURS, name = NULL) +
    scale_shape_manual(values = SEX_SHAPES, guide = "none") +
    labs(
      x = sprintf("PC1 (%.1f%%)", 100 * o$variance[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * o$variance[2]), title = title,
      subtitle = with(o$permanova, paste(
        sprintf("%s R² %.2f, p = %.3f", term[1:3], r2[1:3], p[1:3]),
        collapse = "; "
      ))
    ) +
    theme_figure()
}
r2_page <- ggplot(component_r2, aes(pc, factor(factor, rev(names(codes))), fill = r2)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r2)), size = 2.4) +
  scale_fill_gradient(low = "white", high = "grey25", limits = c(0, 1)) +
  labs(x = NULL, y = NULL, title = "What each component tracks", subtitle = "R² of each factor alone") +
  theme_figure()
plot_index <- write_report(file.path(out$reports, "10_global.pdf"), list(
  "proteome PCA" = pca_page(proteome_scores, proteome, "Proteome PCA"),
  "component R²" = r2_page,
  "phenotype PCA" = pca_page(phenotype_scores, phenotype, "Phenotype PCA, six traits scaled")
))

saveRDS(
  list(
    proteome_scores = proteome_scores, proteome_variance = proteome$variance,
    permanova = proteome$permanova, dispersion_p = proteome$dispersion_p,
    component_r2 = component_r2, mito_vector = mito_vector, technical = samples,
    technical_tests = technical_tests,
    phenotype_scores = phenotype_scores, phenotype_variance = phenotype$variance,
    phenotype_permanova = phenotype$permanova, phenotype_loadings = phenotype_loadings
  ),
  file.path(out$data, "global.rds")
)
write_workbook(
  list(
    proteome_scores = proteome_scores,
    proteome_variance = tibble(component = seq_along(proteome$variance), share = proteome$variance),
    permanova = mutate(proteome$permanova, dispersion_p = c(
      proteome$dispersion_p, rep(NA, nrow(proteome$permanova) - 1)
    )),
    component_r2 = select(component_r2, component, factor, r2),
    mito_vector = mito_vector,
    technical_groups = technical_groups, technical_tests = technical_tests,
    phenotype_scores = phenotype_scores, phenotype_loadings = phenotype_loadings,
    phenotype_permanova = mutate(phenotype$permanova, dispersion_p = c(
      phenotype$dispersion_p, rep(NA, nrow(phenotype$permanova) - 1)
    )),
    plot_index = plot_index
  ),
  c(
    "Proteome PCA scores of the twenty mice, normalised full universe, centred, not scaled.",
    "Share of proteome variance per component.",
    "Proteome PERMANOVA (adonis2, Euclidean, 999 permutations, sequential); betadisper p, row 1.",
    "R² of each design factor, the MitoCarta share and run depth on PC1 to PC4.",
    "MitoCarta share fitted onto PC1-PC2 (envfit): unit direction, R² and permutation p.",
    "Median run depth, protein groups, relative raw signal and MitoCarta share per group.",
    "Each technical variable against treatment, sex and their interaction (lm p).",
    "Phenotype PCA scores, six traits scaled.",
    "Phenotype PCA loadings on PC1 and PC2.",
    "Phenotype PERMANOVA on scaled traits; betadisper p on row 1.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

proteome$permanova
phenotype$permanova

sessionInfo()
