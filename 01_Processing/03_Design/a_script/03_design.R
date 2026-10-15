# The model and the questions asked of it: four group means, one surrogate variable for unwanted
# variation, and five contrasts. 01_Fit reads what this writes and changes neither; the README
# gives the reasons.

source(here::here("00_Setup", "Shared_functions", "helpers.R"))
library(limma)

out <- stage_paths("01_Processing", "03_Design")
inputs <- c(
  proteins = here("01_Processing", "02_Quantification", "c_data", "proteins.rds"),
  filtering = here("01_Processing", "01_Filtering", "c_data", "precursors_filtered.rds")
)
quant <- readRDS(inputs[["proteins"]])
proteins <- quant$proteins
run_qc <- readRDS(inputs[["filtering"]])$run_qc

groups <- model.matrix(~ 0 + group, data = proteins$targets)
colnames(groups) <- sub("^group", "", colnames(groups))
rownames(groups) <- proteins$targets$sample_id

# One surrogate variable (Leek & Storey 2007): the strongest structure left in the normalised full
# proteome once the four group means are removed. One, as the sva vignette fixes for small studies;
# the model then spends one residual df on it. Its sign is set so it rises with run depth.
sv1 <- sva::sva(quant$E_norm$full,
  mod = groups, mod0 = matrix(1, ncol(proteins), 1), n.sv = 1
)$sv[, 1]
run_depth <- log2(run_qc$precursors[match(colnames(proteins), run_qc$sample_id)])
sv1 <- sv1 * sign(cor(sv1, run_depth))
design <- cbind(groups, sv1 = sv1)

contrasts <- makeContrasts(
  PAS_vs_VEH_female = PAS_F - VEH_F,
  PAS_vs_VEH_male = PAS_M - VEH_M,
  treatment_by_sex = (PAS_F - VEH_F) - (PAS_M - VEH_M),
  sex_in_VEH = VEH_M - VEH_F,
  sex_in_PAS = PAS_M - PAS_F,
  levels = design
)
stopifnot(
  identical(colnames(contrasts), CONTRASTS),
  identical(rownames(design), colnames(proteins$E)),
  all.equal(
    contrasts[, "treatment_by_sex"], contrasts[, "sex_in_VEH"] - contrasts[, "sex_in_PAS"]
  )
)

# The contrasts are checked on the group-means model, where each coefficient is a raw group mean.
group_means <- t(apply(proteins$E, 1, \(x) tapply(x, proteins$targets$group, mean)))
manual <- cbind(
  PAS_vs_VEH_female = group_means[, "PAS_F"] - group_means[, "VEH_F"],
  PAS_vs_VEH_male = group_means[, "PAS_M"] - group_means[, "VEH_M"],
  treatment_by_sex = (group_means[, "PAS_F"] - group_means[, "VEH_F"]) -
    (group_means[, "PAS_M"] - group_means[, "VEH_M"]),
  sex_in_VEH = group_means[, "VEH_M"] - group_means[, "VEH_F"],
  sex_in_PAS = group_means[, "PAS_M"] - group_means[, "PAS_F"]
)
fitted <- contrasts.fit(lmFit(proteins$E, groups), contrasts[colnames(groups), ])$coefficients
contrast_check <- tibble(
  contrast = colnames(contrasts),
  coefficients_sum_to_zero = abs(colSums(contrasts)) < 1e-12,
  max_abs_diff_from_group_means = signif(
    apply(abs(fitted - manual[, colnames(contrasts)]), 2, max), 2
  )
)
stopifnot(
  all(contrast_check$coefficients_sum_to_zero),
  contrast_check$max_abs_diff_from_group_means < 1e-10
)

surrogate <- tibble(
  sample_id = colnames(proteins), group = proteins$targets$group, sex = proteins$targets$sex,
  sv1 = sv1, run_depth_log2 = run_depth
)
surrogate_r <- tibble(
  mice = c("all", "female", "male"),
  r_with_run_depth = c(
    cor(sv1, run_depth),
    map_dbl(c("F", "M"), \(s) cor(sv1[surrogate$sex == s], run_depth[surrogate$sex == s]))
  ),
  treatment_p = c(NA, map_dbl(c("F", "M"), \(s) {
    keep <- surrogate$sex == s
    t.test(sv1[keep] ~ proteins$targets$treatment[keep])$p.value
  }))
)
plot_index <- preview_pdf(file.path(out$reports, "03_design.pdf"), list(
  "surrogate variable against run depth" = \() {
    plot(run_depth, sv1,
      pch = c(F = 16, M = 17)[surrogate$sex], col = GROUP_COLOURS[surrogate$group],
      xlab = "Run depth (log2 precursors identified)", ylab = "Surrogate variable 1",
      main = sprintf("r = %.2f", surrogate_r$r_with_run_depth[1])
    )
    legend("topleft", FIGURE_GROUPS,
      col = GROUP_COLOURS[FIGURE_GROUPS], pch = c(16, 16, 17, 17), bty = "n"
    )
  }
))

saveRDS(list(design = design, contrasts = contrasts), file.path(out$data, "design.rds"))

write_workbook(
  list(
    design = as_tibble(design, rownames = "sample_id"),
    contrasts = as_tibble(contrasts, rownames = "group"),
    contrast_check = contrast_check,
    comparisons = tibble(contrast = CONTRASTS, comparison = COMPARISON_OF_CONTRAST[CONTRASTS]),
    groups = count(proteins$targets, group, name = "n"),
    surrogate = surrogate, surrogate_r = surrogate_r, plot_index = plot_index
  ),
  c(
    "Design matrix, one row per mouse: four group means and the surrogate variable sv1.",
    "Contrast matrix: the five questions as group weights.",
    "Each contrast checked against its intended group difference.",
    "Which comparison reads each contrast; the baseline sex contrasts are context only.",
    "Mice per group.",
    "Surrogate variable and run depth per mouse.",
    "Surrogate variable's r with run depth, overall and within sex; Welch p for PAS v vehicle.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

contrast_check

sessionInfo()
