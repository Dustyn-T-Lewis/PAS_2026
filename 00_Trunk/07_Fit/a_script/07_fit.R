#!/usr/bin/env Rscript
# Every model the arcs read, fitted once on all 20 mice. Proteins: limpa's dpcDE with precision
# and sample weights, quantile-normalised per universe. GO:BP, MitoCarta, complex and cluster
# scores and module eigengenes: limma. Each level gets the five contrasts and, for each phenotype,
# a link in females and a link in males: the feature regressed on the z-scored phenotype across the
# ten mice of each sex, PAS and vehicle together (~0 + sex + z:female + z:male). The links ask
# whether features that move with treatment move with the phenotypes too; co-movement, not cause.
# A refit with the pellet's mitochondrial share (centred within group) is kept as a check.
# Phenotype effects are % changes against the same sex's vehicle. Arcs slice these tables.

source(here::here("R", "helpers.R"))
library(limma)
library(limpa)

out <- stage_paths("00_Trunk", "07_Fit")
inputs <- c(
  proteins = here("00_Trunk", "02_Quantification", "c_data", "proteins.rds"),
  design = here("00_Trunk", "03_Design", "c_data", "design.rds"),
  scores = here("00_Trunk", "05_Scores", "c_data", "scores.rds"),
  modules = here("00_Trunk", "06_Modules", "c_data", "modules.rds"),
  phenotypes = input_path("phenotypes.xlsx")
)
quant <- readRDS(inputs[["proteins"]])
d <- readRDS(inputs[["design"]])
sc <- readRDS(inputs[["scores"]])
mods <- readRDS(inputs[["modules"]])
targets <- quant$proteins$targets
stopifnot(identical(rownames(d$design), targets$sample_id))
ph <- read_phenotypes(targets$sample_id)
female <- as.numeric(targets$sex == "F")
pellet_mito <- quant$pct_mito - ave(quant$pct_mito, targets$group)
sex_design <- model.matrix(~ 0 + sex, targets)
colnames(sex_design) <- c("F", "M")
genes <- as_tibble(quant$proteins$genes, rownames = "feature") |>
  transmute(feature, gene = str_split_i(Genes, ";", 1), genes = Genes, is_mito)

# The phenotype design: a mean per sex and the z-scored phenotype's slope within each sex, with
# PAS and vehicle mice together; optionally the pellet's mitochondrial share around its group mean.
link_design <- function(x, adjusted) {
  z <- as.numeric(scale(x))
  m <- cbind(sex_design, link_F = z * female, link_M = z * (1 - female))
  if (adjusted) cbind(m, pellet_mito = pellet_mito) else m
}
link_contrasts <- \(design) makeContrasts(female = link_F, male = link_M, levels = design)

tidy <- function(fit, coefs) {
  map(coefs, \(k) {
    topTable(fit, coef = k, number = Inf, sort.by = "none", confint = TRUE) |>
      as_tibble(rownames = "feature") |>
      transmute(feature, coef = k, estimate = logFC, ci_low = CI.L, ci_high = CI.R,
        AveExpr, t, P.Value, adj.P.Val
      )
  }) |>
    list_rbind()
}

# Proteins, one dpcDE fit per universe and design.
protein_elist <- function(universe) {
  p <- quant$proteins[quant$universes[[universe]], ]
  p$E <- quant$E_norm[[universe]]
  p
}
fit_proteins <- function(design, contrasts) {
  imap(quant$universes, \(keep, universe) {
    dpc <- dpcDE(protein_elist(universe), design, sample.weights = TRUE, plot = FALSE)
    list(dpc = dpc, fit = eBayes(contrasts.fit(dpc, contrasts)))
  })
}
protein_fits <- fit_proteins(d$design, d$contrasts)
protein_contrasts <- imap(protein_fits, \(f, universe) {
  tidy(f$fit, CONTRASTS) |> mutate(universe = universe, .before = 1)
}) |>
  list_rbind() |>
  rename(contrast = coef, logFC = estimate) |>
  left_join(genes, by = "feature") |>
  mutate(pi_score = pi_score(P.Value, logFC))

protein_links <- expand_grid(phenotype = PHENOTYPES, adjusted = c(FALSE, TRUE)) |>
  pmap(\(phenotype, adjusted) {
    design <- link_design(ph[[phenotype]], adjusted)
    imap(fit_proteins(design, link_contrasts(design)), \(f, universe) {
      tidy(f$fit, c("female", "male")) |>
        mutate(universe = universe, phenotype = phenotype, adjusted = adjusted, .before = 1)
    }) |>
      list_rbind()
  }) |>
  list_rbind() |>
  rename(sex = coef) |>
  left_join(genes, by = "feature")

# Scores and eigengenes, plain limma: one value per mouse, no precision weights.
score_sets <- list(
  full = sc$scores$full, mito = sc$scores$mito,
  module = t(as.matrix(mods$eigengenes))
)
level_of <- \(universe, feature) {
  case_when(
    universe == "module" ~ "module",
    sc$collection[feature] == "CORUM" ~ "complex",
    sc$collection[feature] == "STRING" ~ "cluster",
    .default = "pathway"
  )
}
fit_scores <- function(design, contrasts, coefs) {
  imap(score_sets, \(m, universe) {
    stopifnot(identical(colnames(m), targets$sample_id))
    tidy(eBayes(contrasts.fit(lmFit(m, design), contrasts)), coefs) |>
      mutate(
        level = level_of(universe, feature),
        collection = if_else(level == "module", "WGCNA", unname(sc$collection[feature])),
        universe = if_else(universe == "module", "full", universe),
        .before = 1
      )
  }) |>
    list_rbind()
}
score_contrasts <- fit_scores(d$design, d$contrasts, CONTRASTS) |>
  rename(contrast = coef, logFC = estimate)
score_links <- expand_grid(phenotype = PHENOTYPES, adjusted = c(FALSE, TRUE)) |>
  pmap(\(phenotype, adjusted) {
    design <- link_design(ph[[phenotype]], adjusted)
    fit_scores(design, link_contrasts(design), c("female", "male")) |>
      mutate(phenotype = phenotype, adjusted = adjusted, .before = 1)
  }) |>
  list_rbind() |>
  rename(sex = coef)

# Effect sizes beside every test, as Pearson r on the normalised values. For a two-group contrast,
# r with group membership (first group = 1), the r-equivalent of a two-group t; empty for the
# interaction. For a phenotype link, r across the sex's ten mice and, descriptively, within its
# five vehicle and five PAS mice.
group_r <- function(m, contrast) {
  g <- COMPARISONS[[contrast]]
  if (is.null(g)) {
    return(tibble(feature = rownames(m), r_group = NA_real_))
  }
  keep <- targets$group %in% g
  tibble(feature = rownames(m), r_group = cor(t(m[, keep]), as.numeric(targets$group[keep] == g[1]))[, 1])
}
link_r <- function(m, x) {
  map(c(female = "F", male = "M"), \(code) {
    r <- \(keep) cor(t(m[, keep]), x[keep])[, 1]
    in_sex <- targets$sex == code
    tibble(
      feature = rownames(m), r = r(in_sex), r_veh = r(in_sex & targets$treatment == "VEH"),
      r_pas = r(in_sex & targets$treatment == "PAS")
    )
  }) |>
    list_rbind(names_to = "sex")
}
matrices <- c(
  set_names(quant$E_norm, paste0("protein_", names(quant$E_norm))),
  set_names(score_sets, paste0("score_", names(score_sets)))
)
r_groups <- imap(matrices, \(m, name) {
  map(CONTRASTS, \(cn) mutate(group_r(m, cn), contrast = cn)) |>
    list_rbind() |>
    mutate(source = name)
}) |>
  list_rbind()
r_links <- imap(matrices, \(m, name) {
  map(PHENOTYPES, \(ph_name) mutate(link_r(m, ph[[ph_name]]), phenotype = ph_name)) |>
    list_rbind() |>
    mutate(source = name)
}) |>
  list_rbind()
protein_contrasts <- protein_contrasts |>
  left_join(
    filter(r_groups, str_starts(source, "protein")) |> mutate(universe = str_remove(source, "protein_")) |> select(-source),
    by = c("universe", "contrast", "feature")
  )
score_contrasts <- score_contrasts |>
  left_join(filter(r_groups, str_starts(source, "score")) |> select(-source), by = c("contrast", "feature"))
protein_links <- protein_links |>
  left_join(
    filter(r_links, str_starts(source, "protein")) |> mutate(universe = str_remove(source, "protein_")) |> select(-source),
    by = c("universe", "phenotype", "sex", "feature")
  )
score_links <- score_links |>
  left_join(filter(r_links, str_starts(source, "score")) |> select(-source), by = c("phenotype", "sex", "feature"))

# Without weights or covariates the joint model's link equals a regression within one sex.
check <- filter(score_links, !adjusted, universe == "full", sex == "female") |>
  slice_head(n = 1, by = phenotype)
stopifnot(pmap_lgl(check, \(phenotype, feature, estimate, ...) {
  keep <- female == 1
  y <- score_sets$full[feature, keep]
  isTRUE(all.equal(unname(coef(lm(y ~ scale(ph[[phenotype]])[keep]))[2]), estimate))
}))

# Phenotype effects: PAS against the same sex's vehicle as a % change (Welch t on log values,
# exponentiated), and the treatment-by-sex term on the same log scale as female minus male.
phenotype_long <- ph |>
  select(mouse, treatment, sex, all_of(PHENOTYPES)) |>
  pivot_longer(all_of(PHENOTYPES), names_to = "phenotype") |>
  mutate(treatment = factor(treatment, c("VEH", "PAS")), sex = factor(sex, c("F", "M")))
phenotype_effects <- phenotype_long |>
  reframe(
    {
      w <- t.test(log(value[treatment == "PAS"]), log(value[treatment == "VEH"]))
      tibble(
        n_pas = sum(treatment == "PAS"), n_veh = sum(treatment == "VEH"),
        pct_change = 100 * (exp(unname(diff(rev(w$estimate)))) - 1),
        ci_low = 100 * (exp(w$conf.int[1]) - 1), ci_high = 100 * (exp(w$conf.int[2]) - 1),
        P.Value = w$p.value
      )
    },
    .by = c(phenotype, sex)
  )
phenotype_interaction <- phenotype_long |>
  reframe(
    {
      m <- lm(log(value) ~ treatment * sex)
      k <- "treatmentPAS:sexM"
      ci <- confint(m)[k, ]
      tibble(
        pct_female_vs_male = 100 * (exp(-coef(m)[[k]]) - 1),
        ci_low = 100 * (exp(-ci[[2]]) - 1), ci_high = 100 * (exp(-ci[[1]]) - 1),
        P.Value = summary(m)$coefficients[k, "Pr(>|t|)"]
      )
    },
    .by = phenotype
  )

# enrichVolcano study folders the arcs' pathway tests read: results, matrix, weights, samples and
# the contrast expressions rebuilt from the design's contrast matrix.
stopifnot(all(d$contrasts %in% c(-1, 0, 1)))
contrast_expression <- apply(d$contrasts, 2, \(k) {
  paste(c(paste(names(k)[k > 0], collapse = " + "), names(k)[k < 0]), collapse = " - ")
})
iwalk(quant$universes, \(keep, universe) {
  dir <- file.path(out$data, paste0("study_", universe))
  dir.create(dir, showWarnings = FALSE)
  protein_contrasts |>
    filter(universe == !!universe) |>
    transmute(protein = feature, gene, contrast, logFC, t, P.Value, adj.P.Val, AveExpr) |>
    readr::write_csv(file.path(dir, "da_results.csv"))
  as_tibble(quant$E_norm[[universe]], rownames = "protein") |>
    readr::write_csv(file.path(dir, "matrix.csv"))
  as_tibble(protein_fits[[universe]]$dpc$EList$weights, rownames = "protein") |>
    set_names(c("protein", targets$sample_id)) |>
    readr::write_csv(file.path(dir, "weights.csv"))
  tibble(sample = targets$sample_id, group = targets$group) |>
    readr::write_csv(file.path(dir, "samples.csv"))
  tibble(name = names(contrast_expression), expression = unname(contrast_expression)) |>
    readr::write_csv(file.path(dir, "contrasts.csv"))
})

contrast_summary <- bind_rows(
  mutate(protein_contrasts, level = "protein"),
  score_contrasts
) |>
  reframe(
    tier_counts(P.Value, adj.P.Val, pi_score),
    share_true_effects = if (n() >= 100) 1 - propTrueNull(P.Value) else NA_real_,
    .by = c(level, universe, contrast)
  )
link_summary <- bind_rows(mutate(protein_links, level = "protein"), score_links) |>
  reframe(tier_counts(P.Value, adj.P.Val), .by = c(level, universe, phenotype, sex, adjusted))

volcano_page <- protein_contrasts |>
  filter(universe == "full") |>
  mutate(
    contrast = factor(CONTRAST_LABELS[contrast], CONTRAST_LABELS),
    tier = case_when(
      adj.P.Val <= ALPHA ~ "FDR ≤ 0.05", pi_score <= ALPHA ~ "Π ≤ 0.05",
      P.Value <= ALPHA ~ "p ≤ 0.05", .default = "p > 0.05"
    )
  ) |>
  ggplot(aes(logFC, -log10(P.Value), colour = tier)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~contrast, nrow = 1) +
  scale_colour_manual(values = c(
    "p > 0.05" = "grey80", "p ≤ 0.05" = "grey50", "Π ≤ 0.05" = "#E08214",
    "FDR ≤ 0.05" = "black"
  ), name = NULL) +
  labs(
    x = "log2 fold change", y = "-log10 p", title = "Proteins, every contrast",
    subtitle = "Full proteome; dashed line p = 0.05"
  ) +
  theme_figure() +
  theme(legend.position = "bottom")
grid_page <- link_summary |>
  filter(level == "protein", universe == "full") |>
  mutate(
    ratio = p_le_0.05 / expected_by_chance,
    column = factor(
      paste(if_else(sex == "female", "Female", "Male"), if_else(adjusted, "+ mito share", "")),
      c("Female ", "Female + mito share", "Male ", "Male + mito share")
    ),
    trait = factor(PHENOTYPE_LABELS[phenotype], rev(PHENOTYPE_LABELS))
  ) |>
  ggplot(aes(column, trait, fill = ratio)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = str_glue("{p_le_0.05}\nFDR {fdr_le_0.05}")), size = 2.2) +
  scale_fill_gradient(low = "white", high = "grey35", name = "p ≤ 0.05 /\nchance") +
  labs(
    x = NULL, y = NULL, title = "Proteins moving with each phenotype",
    subtitle = "Across PAS and vehicle mice of each sex; 92 at p ≤ 0.05 by chance"
  ) +
  theme_figure()
effects_page <- phenotype_effect_plot(phenotype_effects,
  title = "PAS effect on each phenotype", subtitle = "Welch t on log values within sex, n = 5 v 5"
)
plot_index <- write_report(file.path(out$reports, "07_fit.pdf"), list(
  "protein volcanoes, every contrast" = volcano_page,
  "protein-phenotype link counts" = grid_page,
  "phenotype effects" = effects_page
))

saveRDS(
  list(
    protein_contrasts = protein_contrasts, protein_links = protein_links,
    score_contrasts = score_contrasts, score_links = score_links,
    phenotype_effects = phenotype_effects, phenotype_interaction = phenotype_interaction,
    E_norm = quant$E_norm, scores = score_sets, targets = targets, phenotypes = ph,
    pellet_mito = quant$pct_mito
  ),
  file.path(out$data, "fit.rds")
)
write_workbook(
  list(
    contrast_summary = contrast_summary, link_summary = link_summary,
    protein_contrasts = protein_contrasts,
    protein_links = filter(protein_links, !adjusted),
    score_contrasts = score_contrasts, score_links = filter(score_links, !adjusted),
    phenotype_effects = phenotype_effects, phenotype_interaction = phenotype_interaction,
    plot_index = plot_index
  ),
  c(
    "Per level, universe and contrast: p ≤ 0.05, chance, Π ≤ 0.05, FDR tiers, true-effect share.",
    "Per level, universe, phenotype and sex: p ≤ 0.05, chance, FDR tiers; adjusted = mito check.",
    "Every protein, universe and contrast: log2 FC with 95% CI, moderated t, p, BH FDR, Π, r with group.",
    "Every protein's link with each phenotype per sex: slope (log2 per SD), t, p, FDR, r (pooled, VEH, PAS).",
    "Every GO:BP, MitoCarta, complex, cluster and module score: the five contrasts, r with group.",
    "Every score's link with each phenotype per sex: slope, t, p, FDR, r (pooled, VEH, PAS).",
    "PAS vs vehicle within sex: % change with 95% CI and Welch p.",
    "Treatment by sex on the log scale: female response relative to male, % with 95% CI.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

filter(contrast_summary, level == "protein", universe == "full")

sessionInfo()
