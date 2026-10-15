# Every model the comparisons read, fitted once on all 20 mice. Proteins: limpa's dpcDE with
# precision and sample weights, quantile-normalised per universe. GO:BP, Reactome, MitoCarta,
# complex and cluster scores and module eigengenes: limma. Each level gets the five contrasts and,
# for each phenotype, an association in females and one in males: the feature regressed on the
# phenotype (z-scored across all 20 mice) within the ten mice of each sex, PAS and vehicle together
# (~0 + sex + z:female + z:male). An association is co-movement, not a PAS effect. A refit with the
# pellet's mitochondrial share (centred within group) is kept as a check. Phenotype effects are %
# changes against the same sex's vehicle. Comparisons slice these tables.

source(here::here("00_Setup", "Shared_functions", "helpers.R"))
library(limma)
library(limpa)

out <- stage_paths("02_Statistics", "01_Fit")
inputs <- c(
  proteins = here("01_Processing", "02_Quantification", "c_data", "proteins.rds"),
  design = here("01_Processing", "03_Design", "c_data", "design.rds"),
  scores = here("01_Processing", "05_Scores", "c_data", "scores.rds"),
  modules = here("01_Processing", "06_Modules", "c_data", "modules.rds"),
  phenotypes = input_path("phenotypes.xlsx")
)
quant <- readRDS(inputs[["proteins"]])
d <- readRDS(inputs[["design"]])
sc <- readRDS(inputs[["scores"]])
mods <- readRDS(inputs[["modules"]])
targets <- quant$proteins$targets
stopifnot(identical(rownames(d$design), targets$sample_id))
ph <- readxl::read_excel(inputs[["phenotypes"]], "phenotypes")
ph <- ph[match(targets$sample_id, ph$mouse), ]
stopifnot(identical(ph$mouse, targets$sample_id), !anyNA(ph[PHENOTYPES]), all(ph[PHENOTYPES] > 0))
female <- as.numeric(targets$sex == "F")
pellet_mito <- quant$pct_mito - ave(quant$pct_mito, targets$group)
sex_design <- model.matrix(~ 0 + sex, targets)
colnames(sex_design) <- c("F", "M")
genes <- as_tibble(quant$proteins$genes, rownames = "feature") |>
  transmute(feature, gene = first_symbol(Genes), genes = Genes, is_mito)

# The phenotype design: a mean per sex and the z-scored phenotype's slope within each sex, with
# PAS and vehicle mice together; optionally the pellet's mitochondrial share around its group mean.
association_design <- function(x, adjusted) {
  z <- as.numeric(scale(x))
  m <- cbind(sex_design, slope_F = z * female, slope_M = z * (1 - female))
  if (adjusted) cbind(m, pellet_mito = pellet_mito) else m
}
slope_contrasts <- \(design) makeContrasts(female = slope_F, male = slope_M, levels = design)

tidy <- function(fit, coefs) {
  map(coefs, \(k) {
    topTable(fit, coef = k, number = Inf, sort.by = "none", confint = TRUE) |>
      as_tibble(rownames = "feature") |>
      transmute(feature,
        coef = k, estimate = logFC, ci_low = CI.L, ci_high = CI.R,
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
  map(set_names(names(quant$universes)), \(universe) {
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

protein_associations <- expand_grid(phenotype = PHENOTYPES, adjusted = c(FALSE, TRUE)) |>
  pmap(\(phenotype, adjusted) {
    design <- association_design(ph[[phenotype]], adjusted)
    imap(fit_proteins(design, slope_contrasts(design)), \(f, universe) {
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
score_associations <- expand_grid(phenotype = PHENOTYPES, adjusted = c(FALSE, TRUE)) |>
  pmap(\(phenotype, adjusted) {
    design <- association_design(ph[[phenotype]], adjusted)
    fit_scores(design, slope_contrasts(design), c("female", "male")) |>
      mutate(phenotype = phenotype, adjusted = adjusted, .before = 1)
  }) |>
  list_rbind() |>
  rename(sex = coef)

# Effect sizes beside every test, as Pearson r on the normalised values. For a two-group contrast,
# r with group membership (first group = 1), the r-equivalent of a two-group t; empty for the
# interaction. For a phenotype association, r across the sex's ten mice and, descriptively, within
# its five vehicle and five PAS mice.
GROUP_PAIRS <- list(
  PAS_vs_VEH_female = c("PAS_F", "VEH_F"),
  PAS_vs_VEH_male = c("PAS_M", "VEH_M"),
  sex_in_VEH = c("VEH_M", "VEH_F"),
  sex_in_PAS = c("PAS_M", "PAS_F")
)
group_r <- function(m, contrast) {
  g <- GROUP_PAIRS[[contrast]]
  if (is.null(g)) {
    return(tibble(feature = rownames(m), r_group = NA_real_))
  }
  keep <- targets$group %in% g
  first <- as.numeric(targets$group[keep] == g[1])
  tibble(feature = rownames(m), r_group = cor(t(m[, keep]), first)[, 1])
}
association_r <- function(m, x) {
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
r_associations <- imap(matrices, \(m, name) {
  map(PHENOTYPES, \(ph_name) mutate(association_r(m, ph[[ph_name]]), phenotype = ph_name)) |>
    list_rbind() |>
    mutate(source = name)
}) |>
  list_rbind()
# Score, module and protein names are unique within their source, so each contrast row takes
# exactly one r; an association's r is shared by its unadjusted row and its mito-share check.
protein_r <- \(r) {
  filter(r, str_starts(source, "protein")) |>
    mutate(universe = str_remove(source, "protein_")) |>
    select(-source)
}
score_r <- \(r) select(filter(r, str_starts(source, "score")), -source)
join_r <- \(x, y, by, relationship = "one-to-one") {
  left_join(x, y, by = by, relationship = relationship, unmatched = "error")
}
protein_contrasts <- protein_contrasts |>
  join_r(protein_r(r_groups), c("universe", "contrast", "feature"))
score_contrasts <- join_r(score_contrasts, score_r(r_groups), c("contrast", "feature"))
protein_associations <- protein_associations |>
  join_r(protein_r(r_associations), c("universe", "phenotype", "sex", "feature"), "many-to-one")
score_associations <- score_associations |>
  join_r(score_r(r_associations), c("phenotype", "sex", "feature"), "many-to-one")

# Without weights or covariates the joint model's association equals a regression within one sex.
check <- filter(score_associations, !adjusted, universe == "full", sex == "female") |>
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

# enrichVolcano study folders 02_Enrichment reads: results, matrix, weights, samples and the
# contrast expressions rebuilt from the design's contrast matrix.
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
association_summary <- protein_associations |>
  mutate(level = "protein") |>
  bind_rows(score_associations) |>
  reframe(tier_counts(P.Value, adj.P.Val), .by = c(level, universe, phenotype, sex, adjusted))

volcano_page <- protein_contrasts |>
  filter(universe == "full") |>
  mutate(
    contrast = factor(CONTRAST_LABELS[contrast], CONTRAST_LABELS),
    tier = case_when(
      adj.P.Val <= ALPHA ~ "FDR ≤ 0.05", pi_score <= ALPHA ~ "Π ≤ 0.05",
      P.Value <= ALPHA ~ "p ≤ 0.05",
      .default = "p > 0.05"
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
grid_page <- association_summary |>
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
    x = NULL, y = NULL, title = "Proteins associated with each phenotype",
    subtitle = str_glue(
      "Across PAS and vehicle mice of each sex; ",
      "{round(ALPHA * sum(quant$universes$full))} at p ≤ 0.05 by chance"
    )
  ) +
  theme_figure()
effects_page <- phenotype_effect_plot(phenotype_effects,
  title = "Phenotype effects of PAS",
  subtitle = str_glue(
    "Welch t on log values within sex, n = {min(phenotype_effects$n_pas)} v ",
    "{min(phenotype_effects$n_veh)}"
  )
)
plot_index <- write_report(file.path(out$reports, "01_fit.pdf"), list(
  "protein volcanoes, every contrast" = volcano_page,
  "protein-phenotype association counts" = grid_page,
  "phenotype effects" = effects_page
))

saveRDS(
  list(
    protein_contrasts = protein_contrasts, protein_associations = protein_associations,
    score_contrasts = score_contrasts, score_associations = score_associations,
    phenotype_effects = phenotype_effects, phenotype_interaction = phenotype_interaction,
    E_norm = quant$E_norm, scores = score_sets, targets = targets, phenotypes = ph,
    pct_mito = quant$pct_mito
  ),
  file.path(out$data, "fit.rds")
)
write_workbook(
  list(
    contrast_summary = contrast_summary, association_summary = association_summary,
    protein_contrasts = protein_contrasts,
    protein_associations = filter(protein_associations, !adjusted),
    score_contrasts = score_contrasts, score_associations = filter(score_associations, !adjusted),
    phenotype_effects = phenotype_effects, phenotype_interaction = phenotype_interaction,
    plot_index = plot_index
  ),
  c(
    "Per level, universe and contrast: p ≤ 0.05, chance, Π ≤ 0.05, FDR tiers, true-effect share.",
    "Per level, universe, phenotype and sex: p ≤ 0.05, chance, FDR tiers; adjusted = mito check.",
    paste(
      "Every protein, universe and contrast: log2 FC with 95% CI, moderated t, p, BH FDR, Π,",
      "r with group."
    ),
    paste(
      "Every protein's association with each phenotype per sex: slope (log2 per SD across all",
      "20 mice), t, p, FDR, r (pooled, VEH, PAS)."
    ),
    paste(
      "Every GO:BP, Reactome, MitoCarta, complex, cluster and module score: the five contrasts,",
      "r with group."
    ),
    paste(
      "Every score's association with each phenotype per sex: slope, t, p, FDR,",
      "r (pooled, VEH, PAS)."
    ),
    "Phenotype effects: PAS vs vehicle within sex, % change with 95% CI and Welch p.",
    paste(
      "Interaction (treatment_by_sex) on the log scale: female response relative to male,",
      "% with 95% CI."
    ),
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

filter(contrast_summary, level == "protein", universe == "full")

sessionInfo()
