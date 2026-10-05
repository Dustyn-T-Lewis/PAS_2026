#!/usr/bin/env Rscript
# What shapes the whole proteome: PCA of the twenty mice on the normalised full universe (centred,
# not scaled), PERMANOVA (vegan::adonis2, Euclidean, 999 permutations, terms added in sequence) for
# sex, treatment and their interaction, a dispersion check (betadisper by group), and how much of
# each of the first four components sex, treatment, the interaction and the MitoCarta share explain
# (R², each on its own; the design codes are orthogonal in this balanced 2 x 2). The MitoCarta share
# is fitted onto PC1-PC2 as a vector (vegan::envfit, 999 permutations).

source(here::here("R", "helpers.R"))

N_COMPONENTS <- 4
PERMUTATIONS <- 999

out <- stage_paths("05_Summary", "05_Global")
inputs <- c(
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
de <- readRDS(inputs[["contrasts"]])
samples <- as_tibble(de$targets) |>
  select(sample_id, treatment, sex, group) |>
  mutate(pct_mito = de$pct_mito[sample_id])
x <- t(de$E_norm$full[, samples$sample_id])

pca <- prcomp(x, center = TRUE, scale. = FALSE)
variance <- pca$sdev^2 / sum(pca$sdev^2)
scores <- as_tibble(pca$x[, seq_len(N_COMPONENTS)]) |>
  bind_cols(samples, .name_repair = "minimal")

distances <- dist(x)
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
dispersion_p <- dispersion$tab[1, "Pr(>F)"]

set.seed(SEED)
mito_fit <- vegan::envfit(pca$x[, 1:2], samples["pct_mito"], permutations = PERMUTATIONS)
mito_vector <- tibble(
  PC1 = mito_fit$vectors$arrows[1, 1], PC2 = mito_fit$vectors$arrows[1, 2],
  r2 = unname(mito_fit$vectors$r), p = unname(mito_fit$vectors$pvals)
)

codes <- transmute(samples,
  Sex = if_else(sex == "M", 1, -1),
  Treatment = if_else(treatment == "PAS", 1, -1),
  `Sex × treatment` = Sex * Treatment,
  `MitoCarta share` = pct_mito
)
component_r2 <- expand_grid(component = seq_len(N_COMPONENTS), factor = names(codes)) |>
  mutate(
    r2 = map2_dbl(component, factor, \(k, f) cor(pca$x[, k], codes[[f]])^2),
    pc = factor(sprintf("PC%d\n(%.1f%%)", component, 100 * variance[component]))
  )

saveRDS(
  list(
    scores = scores, variance = variance, permanova = permanova,
    dispersion_p = dispersion_p, component_r2 = component_r2, mito_vector = mito_vector
  ),
  file.path(out$data, "global.rds")
)
write_workbook(
  list(
    scores = scores,
    variance = tibble(component = seq_along(variance), share = variance),
    permanova = mutate(permanova, dispersion_p = c(dispersion_p, rep(NA, nrow(permanova) - 1))),
    component_r2 = select(component_r2, component, factor, r2),
    mito_vector = mito_vector
  ),
  c(
    "PCA scores of the twenty mice, normalised full universe, centred, not scaled.",
    "Share of variance per component.",
    "PERMANOVA (adonis2, Euclidean, 999 permutations, sequential terms); betadisper p on row 1.",
    "R² of each design factor and the MitoCarta share on PC1 to PC4.",
    "MitoCarta share fitted onto PC1-PC2 (envfit): unit direction, R² and permutation p."
  ),
  out$workbook, inputs
)

permanova
dispersion_p
sessionInfo()
