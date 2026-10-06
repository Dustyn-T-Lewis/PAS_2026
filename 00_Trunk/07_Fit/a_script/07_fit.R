#!/usr/bin/env Rscript
# Every model the arcs read, fitted once on all 20 mice. Proteins: limpa's dpcDE with precision
# and sample weights, quantile-normalised per universe. Pathway, complex and cluster scores and
# module eigengenes: limma. Each level gets the five contrasts and, for each phenotype, a female
# slope, a male slope and their difference from one model (~0 + group + z:female + z:male +
# within-group pellet mitochondrial share). Phenotype effects are % changes against the same sex's
# vehicle. Arcs slice these tables; nothing downstream refits them.

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
genes <- as_tibble(quant$proteins$genes, rownames = "feature") |>
  transmute(feature, gene = str_split_i(Genes, ";", 1), genes = Genes, is_mito)

# The phenotype design: one z-scored phenotype, its slope within females and within males beside
# the four group means, and optionally the pellet's mitochondrial share around its group mean.
slope_design <- function(x, adjusted) {
  z <- as.numeric(scale(x))
  m <- cbind(d$design, slope_F = z * female, slope_M = z * (1 - female))
  if (adjusted) cbind(m, pellet_mito = pellet_mito) else m
}
slope_contrasts <- \(design) {
  makeContrasts(
    female = slope_F, male = slope_M, difference = slope_M - slope_F,
    levels = design
  )
}

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

protein_slopes <- expand_grid(phenotype = PHENOTYPES, adjusted = c(TRUE, FALSE)) |>
  pmap(\(phenotype, adjusted) {
    design <- slope_design(ph[[phenotype]], adjusted)
    imap(fit_proteins(design, slope_contrasts(design)), \(f, universe) {
      tidy(f$fit, c("female", "male", "difference")) |>
        mutate(universe = universe, phenotype = phenotype, adjusted = adjusted, .before = 1)
    }) |>
      list_rbind()
  }) |>
  list_rbind() |>
  rename(slope = coef) |>
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
        universe = if_else(universe == "module", "full", universe),
        level = level_of(universe, feature),
        collection = if_else(level == "module", "WGCNA", unname(sc$collection[feature])),
        .before = 1
      )
  }) |>
    list_rbind()
}
score_contrasts <- fit_scores(d$design, d$contrasts, CONTRASTS) |>
  rename(contrast = coef, logFC = estimate)
score_slopes <- expand_grid(phenotype = PHENOTYPES, adjusted = c(TRUE, FALSE)) |>
  pmap(\(phenotype, adjusted) {
    design <- slope_design(ph[[phenotype]], adjusted)
    fit_scores(design, slope_contrasts(design), c("female", "male", "difference")) |>
      mutate(phenotype = phenotype, adjusted = adjusted, .before = 1)
  }) |>
  list_rbind() |>
  rename(slope = coef)

# The female slope plus the difference is the male slope, for every feature and model.
check <- bind_rows(protein_slopes, score_slopes) |>
  select(any_of(c("universe", "level")), phenotype, adjusted, feature, slope, estimate) |>
  pivot_wider(names_from = slope, values_from = estimate)
stopifnot(max(abs(check$female + check$difference - check$male)) < 1e-8)

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
slope_summary <- bind_rows(mutate(protein_slopes, level = "protein"), score_slopes) |>
  reframe(tier_counts(P.Value, adj.P.Val), .by = c(level, universe, phenotype, slope, adjusted))

volcano_page <- protein_contrasts |>
  filter(universe == "full") |>
  mutate(
    contrast = factor(CONTRAST_LABELS[contrast], CONTRAST_LABELS),
    tier = case_when(
      adj.P.Val <= ALPHA ~ "FDR ≤ 0.05", P.Value <= ALPHA & pi_score <= ALPHA ~ "p and Π ≤ 0.05",
      P.Value <= ALPHA ~ "p ≤ 0.05", .default = "p > 0.05"
    )
  ) |>
  ggplot(aes(logFC, -log10(P.Value), colour = tier)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~contrast, nrow = 1) +
  scale_colour_manual(values = c(
    "p > 0.05" = "grey80", "p ≤ 0.05" = "grey50", "p and Π ≤ 0.05" = "#E08214",
    "FDR ≤ 0.05" = "black"
  ), name = NULL) +
  labs(
    x = "log2 fold change", y = "-log10 p", title = "Proteins, every contrast",
    subtitle = "Full proteome; dashed line p = 0.05"
  ) +
  theme_figure() +
  theme(legend.position = "bottom")
grid_page <- slope_summary |>
  filter(level == "protein", universe == "full", adjusted) |>
  mutate(ratio = p_le_0.05 / expected_by_chance, slope = factor(slope, c("female", "male", "difference"))) |>
  ggplot(aes(slope, phenotype, fill = ratio)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = str_glue("{p_le_0.05}\nFDR {fdr_le_0.05}")), size = 2.2) +
  scale_fill_gradient(low = "white", high = "grey35", name = "p ≤ 0.05 /\nchance") +
  labs(
    x = NULL, y = NULL, title = "Proteins tracking each phenotype",
    subtitle = "Per-sex slopes, adjusted for pellet mitochondrial share; 92 at p ≤ 0.05 by chance"
  ) +
  theme_figure()
effects_page <- phenotype_effects |>
  mutate(sex = if_else(sex == "F", "Female", "Male")) |>
  ggplot(aes(pct_change, phenotype, colour = sex)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high),
    width = 0, orientation = "y",
    position = position_dodge(width = 0.5)
  ) +
  geom_point(aes(shape = P.Value <= ALPHA), position = position_dodge(width = 0.5), size = 1.8) +
  scale_colour_manual(values = SEX_COLOURS, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c("p > 0.05", "p ≤ 0.05"), name = NULL) +
  labs(
    x = "PAS vs vehicle, % change (95% CI)", y = NULL, title = "PAS effect on each phenotype",
    subtitle = "Welch t on log values within sex, n = 5 v 5"
  ) +
  theme_figure()
plot_index <- write_report(file.path(out$reports, "07_fit.pdf"), list(
  "protein volcanoes, every contrast" = volcano_page,
  "protein-phenotype slope counts" = grid_page,
  "phenotype effects" = effects_page
))

saveRDS(
  list(
    protein_contrasts = protein_contrasts, protein_slopes = protein_slopes,
    score_contrasts = score_contrasts, score_slopes = score_slopes,
    phenotype_effects = phenotype_effects, phenotype_interaction = phenotype_interaction,
    E_norm = quant$E_norm, scores = score_sets, targets = targets, phenotypes = ph,
    pellet_mito = quant$pct_mito
  ),
  file.path(out$data, "fit.rds")
)
write_workbook(
  list(
    contrast_summary = contrast_summary, slope_summary = slope_summary,
    protein_contrasts = protein_contrasts,
    protein_slopes = filter(protein_slopes, adjusted),
    score_contrasts = score_contrasts, score_slopes = filter(score_slopes, adjusted),
    phenotype_effects = phenotype_effects, phenotype_interaction = phenotype_interaction,
    plot_index = plot_index
  ),
  c(
    "Per level, universe and contrast: p ≤ 0.05, chance, p and Π ≤ 0.05, FDR tiers, true-effect share.",
    "Per level, universe, phenotype, slope and adjustment: p ≤ 0.05, chance, FDR tiers.",
    "Every protein, universe and contrast: log2 FC with 95% CI, moderated t, p, BH FDR, Π.",
    "Every protein's female slope, male slope and difference per phenotype (log2 per SD), adjusted.",
    "Every pathway, complex, cluster and module score: the five contrasts.",
    "Every score's female slope, male slope and difference per phenotype, adjusted.",
    "PAS vs vehicle within sex: % change with 95% CI and Welch p.",
    "Treatment by sex on the log scale: female response relative to male, % with 95% CI.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

filter(contrast_summary, level == "protein", universe == "full")

sessionInfo()
