# Every model the comparisons read, fitted once on all mice on 01_Processing/03_Design's design:
# proteins by limpa's dpcDE per universe, scores by limma, eigengenes by ordinary t, each with the
# five contrasts and a per-sex association with each phenotype. Two association models ask two
# questions: pooled, whether a feature co-varies with the phenotype across PAS and vehicle mice;
# treatment, whether it tracks the phenotype among mice given the same treatment.

source(here::here("00_Setup", "Shared_functions", "helpers.R"))
library(limpa)

out <- stage_paths("02_Statistics", "01_Fit")
inputs <- c(
  proteins = output_path("01_Processing", "02_Quantification", "proteins.rds"),
  design = output_path("01_Processing", "03_Design", "design.rds"),
  scores = output_path("01_Processing", "05_Scores", "scores.rds"),
  modules = output_path("01_Processing", "06_Modules", "modules.rds"),
  phenotypes = input_path("phenotypes.xlsx")
)
quant <- readRDS(inputs[["proteins"]])
d <- readRDS(inputs[["design"]])
sc <- readRDS(inputs[["scores"]])
mods <- readRDS(inputs[["modules"]])
targets <- quant$proteins$targets
ph <- readxl::read_excel(inputs[["phenotypes"]], "phenotypes")
ph <- ph[match(targets$sample_id, ph$mouse), ]
female <- as.numeric(targets$sex == "F")
pellet_mito <- quant$pct_mito - ave(quant$pct_mito, targets$group)
sex_design <- model.matrix(~ 0 + sex, targets)
colnames(sex_design) <- c("F", "M")
group_design <- model.matrix(~ 0 + group, targets)
sv1 <- d$design[, "sv1"]
ASSOCIATION_MODELS <- c("pooled", "treatment", "pooled_mito")
# propTrueNull reads the share of true nulls off the p-value distribution, which a few dozen tests
# cannot shape.
PI0_MIN_TESTS <- 100
genes <- as_tibble(quant$proteins$genes, rownames = "feature") |>
  transmute(feature, gene = first_symbol(Genes), genes = Genes, is_mito)

# The phenotype design: a mean per sex (pooled) or per group (treatment) and the z-scored
# phenotype's slope within each sex; the check adds the pellet's mitochondrial share around its
# group mean to the pooled model.
association_design <- function(x, model) {
  z <- as.numeric(scale(x))
  means <- if (model == "treatment") group_design else sex_design
  m <- cbind(means,
    slope_F = z * female, slope_M = z * (1 - female),
    sv1_F = sv1 * female, sv1_M = sv1 * (1 - female)
  )
  if (model == "pooled_mito") cbind(m, pellet_mito = pellet_mito) else m
}
slope_contrasts <- \(design) makeContrasts(female = slope_F, male = slope_M, levels = design)

tidy <- function(fit, coefs) {
  map(coefs, \(k) {
    topTable(fit, coef = k, number = Inf, sort.by = "none", confint = TRUE) |>
      as_tibble(rownames = "feature") |>
      transmute(feature,
        coef = k, estimate = logFC, ci_low = CI.L, ci_high = CI.R,
        AveExpr, t, df = fit$df.total, r_equiv = t / sqrt(t^2 + df), P.Value, adj.P.Val
      )
  }) |>
    list_rbind()
}

# Proteins: one dpcDE fit per universe, design and contrast. contrasts.fit() takes the correlation
# between coefficients as common to all proteins, which is exact only for an orthogonal design or
# equal weights (limma's help page); limpa's weights differ by protein and the surrogate variable
# makes the design non-orthogonal. So each contrast is refitted as the first coefficient of an
# equivalent design (limma::contrastAsCoef): the fitted values and weights are unchanged and the
# contrast's standard error is exact.
protein_elist <- function(universe) {
  p <- quant$proteins[quant$universes[[universe]], ]
  p$E <- quant$E_norm[[universe]]
  p
}
fit_proteins <- function(design, contrasts) {
  map(set_names(names(quant$universes)), \(universe) {
    y <- protein_elist(universe)
    map(set_names(colnames(contrasts)), \(cn) {
      reparametrised <- contrastAsCoef(design, contrasts[, cn])$design
      eBayes(dpcDE(y, reparametrised, sample.weights = TRUE, plot = FALSE))
    })
  })
}
tidy_proteins <- function(fits) {
  imap(fits, \(by_contrast, universe) {
    imap(by_contrast, \(f, cn) mutate(tidy(f, 1), coef = cn)) |>
      list_rbind() |>
      mutate(universe = universe, .before = 1)
  }) |>
    list_rbind()
}
protein_fits <- fit_proteins(d$design, d$contrasts)
# Two calibrated companions to Π, kept as columns: TREAT tests whether a change exceeds 20%
# (McCarthy & Smyth 2009); ashr shrinks each fold change by its standard error and gives the
# local false sign rate, the probability that the direction is wrong (Stephens 2017).
TREAT_LFC <- log2(1.2)
treat_tables <- imap(protein_fits, \(by_contrast, universe) {
  imap(by_contrast, \(f, cn) {
    topTreat(treat(f, lfc = TREAT_LFC), coef = 1, number = Inf, sort.by = "none") |>
      as_tibble(rownames = "feature") |>
      transmute(
        universe = universe, feature, contrast = cn, treat_p = P.Value, treat_fdr = adj.P.Val
      )
  }) |>
    list_rbind()
}) |>
  list_rbind()
shrink <- function(x) {
  fit <- ashr::ash(x$logFC, x$logFC / x$t, df = median(x$df), mixcompdist = "halfuniform")
  mutate(x, logFC_shrunk = ashr::get_pm(fit), lfsr = ashr::get_lfsr(fit))
}
protein_contrasts <- tidy_proteins(protein_fits) |>
  rename(contrast = coef, logFC = estimate) |>
  left_join(genes, by = "feature") |>
  mutate(pi_score = pi_score(P.Value, logFC)) |>
  left_join(treat_tables, by = c("universe", "feature", "contrast")) |>
  group_by(universe, contrast) |>
  group_modify(\(x, key) shrink(x)) |>
  ungroup() |>
  relocate(universe, feature, contrast)

# Does PAS change a protein in either sex? The moderated F over the two PAS contrasts, both made
# the first coefficients of an equivalent design so the test stays exact under limpa's weights.
# It covers the interaction too, which is their difference. The response pattern reads the signs.
pas_contrasts <- d$contrasts[, c("PAS_vs_VEH_female", "PAS_vs_VEH_male")]
protein_ftest <- map(set_names(names(quant$universes)), \(universe) {
  reparametrised <- contrastAsCoef(d$design, pas_contrasts)$design
  f <- eBayes(dpcDE(protein_elist(universe), reparametrised, sample.weights = TRUE, plot = FALSE))
  topTable(f, coef = 1:2, number = Inf, sort.by = "none") |>
    as_tibble(rownames = "feature") |>
    select(feature, F_stat = "F", P.Value, adj.P.Val) |>
    mutate(universe = universe, .before = 1)
}) |>
  list_rbind() |>
  left_join(
    protein_contrasts |>
      filter(contrast %in% colnames(pas_contrasts)) |>
      select(universe, feature, contrast, logFC) |>
      pivot_wider(names_from = contrast, values_from = logFC, names_prefix = "logFC_"),
    by = c("universe", "feature")
  ) |>
  mutate(pattern = case_when(
    logFC_PAS_vs_VEH_female > 0 & logFC_PAS_vs_VEH_male > 0 ~ "up in both",
    logFC_PAS_vs_VEH_female < 0 & logFC_PAS_vs_VEH_male < 0 ~ "down in both",
    logFC_PAS_vs_VEH_female > 0 ~ "female up, male down",
    .default = "female down, male up"
  )) |>
  left_join(genes, by = "feature")

association_fits <- expand_grid(phenotype = PHENOTYPES, model = ASSOCIATION_MODELS) |>
  mutate(
    design = map2(phenotype, model, \(x, m) association_design(ph[[x]], m)),
    fits = map(design, \(x) fit_proteins(x, slope_contrasts(x)))
  )
protein_associations <- association_fits |>
  pmap(\(phenotype, model, design, fits) {
    mutate(tidy_proteins(fits), phenotype = phenotype, model = model, .before = 1)
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
# Ordinary least-squares t on the residual df, for the module eigengenes: eBayes estimates its prior
# from the spread of variances across features, and eight standardised eigengenes cannot inform
# that; it put the prior df at its cap and gave every module about 120 df.
tidy_ordinary <- function(fit, coefs) {
  map(coefs, \(k) {
    estimate <- unname(fit$coefficients[, k])
    se <- unname(fit$stdev.unscaled[, k] * fit$sigma)
    t <- estimate / se
    margin <- qt(0.975, fit$df.residual) * se
    p <- 2 * pt(-abs(t), fit$df.residual)
    tibble(
      feature = rownames(fit$coefficients), coef = k, estimate = estimate,
      ci_low = estimate - margin, ci_high = estimate + margin, AveExpr = unname(fit$Amean), t = t,
      df = fit$df.residual, r_equiv = t / sqrt(t^2 + df), P.Value = p,
      adj.P.Val = p.adjust(p, "BH")
    )
  }) |>
    list_rbind()
}
fit_scores <- function(design, contrasts, coefs) {
  imap(score_sets, \(m, universe) {
    f <- contrasts.fit(lmFit(m, design), contrasts)
    tidied <- if (universe == "module") tidy_ordinary(f, coefs) else tidy(eBayes(f), coefs)
    tidied |>
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
score_associations <- expand_grid(phenotype = PHENOTYPES, model = ASSOCIATION_MODELS) |>
  pmap(\(phenotype, model) {
    design <- association_design(ph[[phenotype]], model)
    fit_scores(design, slope_contrasts(design), c("female", "male")) |>
      mutate(phenotype = phenotype, model = model, .before = 1)
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
# exactly one r; an association's r is shared by its rows under each model.
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

contrast_summary <- bind_rows(
  mutate(protein_contrasts, level = "protein"),
  score_contrasts
) |>
  reframe(
    tier_counts(P.Value, adj.P.Val, pi_score),
    share_true_effects = if (n() >= PI0_MIN_TESTS) 1 - propTrueNull(P.Value) else NA_real_,
    .by = c(level, universe, contrast)
  )
association_summary <- protein_associations |>
  mutate(level = "protein") |>
  bind_rows(score_associations) |>
  reframe(tier_counts(P.Value, adj.P.Val), .by = c(level, universe, phenotype, sex, model))

sa_page <- \(universe) \() plotSA(protein_fits[[universe]][[1]], main = universe)
volcano <- \(cn) {
  \() {
    f <- protein_fits$full[[cn]]
    fdr <- filter(protein_contrasts, universe == "full", contrast == cn)$adj.P.Val
    volcanoplot(f,
      coef = 1, highlight = sum(fdr <= ALPHA),
      names = first_symbol(f$genes$Genes), main = CONTRAST_LABELS[[cn]]
    )
  }
}
plot_index <- preview_pdf(file.path(out$reports, "01_fit.pdf"), c(
  list(
    "protein variance against abundance, full proteome" = sa_page("full"),
    "protein variance against abundance, MitoCarta" = sa_page("mito")
  ),
  set_names(
    map(COMPARISON_CONTRASTS, volcano), paste("volcano,", CONTRAST_LABELS[COMPARISON_CONTRASTS])
  )
))

saveRDS(
  list(
    protein_contrasts = protein_contrasts, protein_ftest = protein_ftest,
    protein_associations = protein_associations, score_contrasts = score_contrasts,
    score_associations = score_associations, association_summary = association_summary,
    phenotype_effects = phenotype_effects, phenotype_interaction = phenotype_interaction,
    E_norm = quant$E_norm, scores = score_sets, targets = targets, phenotypes = ph,
    pct_mito = quant$pct_mito,
    association_designs = select(filter(association_fits, model != "pooled_mito"), -fits)
  ),
  file.path(out$data, "fit.rds")
)
write_workbook(
  list(
    contrast_summary = contrast_summary, association_summary = association_summary,
    protein_contrasts = protein_contrasts, protein_ftest = protein_ftest,
    protein_associations = filter(protein_associations, model != "pooled_mito"),
    score_contrasts = score_contrasts,
    score_associations = filter(score_associations, model != "pooled_mito"),
    phenotype_effects = phenotype_effects, phenotype_interaction = phenotype_interaction,
    plot_index = plot_index
  ),
  c(
    "Per level, universe and contrast: p ≤ 0.05, chance, Π ≤ 0.05, FDR tiers, true-effect share.",
    paste(
      "Per level, universe, phenotype, sex and model: p ≤ 0.05, chance, FDR tiers;",
      "pooled_mito is the pellet mito-share check."
    ),
    paste(
      "Every protein, universe and contrast: log2 FC with 95% CI, moderated t, df, p, BH FDR, Π,",
      "TREAT p and FDR for a change above 20%, ashr shrunken log2 FC and local false sign rate,",
      "r with group."
    ),
    paste(
      "Every protein: moderated F over PAS vs VEH in females and in males (PAS in either sex),",
      "p, BH FDR, both log2 FC and the response pattern from their signs."
    ),
    paste(
      "Every protein's association with each phenotype per sex and model (pooled, treatment):",
      "slope (log2 per SD across all 20 mice), t, p, FDR, r (pooled, VEH, PAS)."
    ),
    paste(
      "Every pathway (GO:BP, Reactome, Hallmark, KEGG, MitoCarta), complex, cluster and module",
      "score: the five contrasts, r with group."
    ),
    paste(
      "Every score's association with each phenotype per sex and model: slope, t, p, FDR,",
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
tibble(
  phenotype_rows_matched = sum(!is.na(ph$mouse)), phenotype_na = sum(is.na(ph[PHENOTYPES])),
  phenotype_le_0 = sum(ph[PHENOTYPES] <= 0, na.rm = TRUE)
)
