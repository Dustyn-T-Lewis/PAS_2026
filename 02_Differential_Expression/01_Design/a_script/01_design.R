#!/usr/bin/env Rscript
# The model and the questions asked of it. 02_Differential reads what this writes and changes
# neither; the README gives the reasons.

library(here)
library(limma)
library(dplyr)
library(tibble)
source(here::here("R", "helpers.R"))

QUANT_DATA <- here("01_Preprocess", "02_Quantification", "c_data")
DESIGN_DATA <- here("02_Differential_Expression", "01_Design", "c_data")

quant <- readRDS(file.path(QUANT_DATA, "proteins.rds"))
proteins <- quant$proteins

design <- model.matrix(~ 0 + group, data = proteins$targets)
colnames(design) <- sub("^group", "", colnames(design))
rownames(design) <- proteins$targets$sample_id

contrasts <- makeContrasts(
  PAS_vs_VEH_female = PAS_F - VEH_F,
  PAS_vs_VEH_male = PAS_M - VEH_M,
  treatment_by_sex = (PAS_F - VEH_F) - (PAS_M - VEH_M),
  sex_in_VEH = VEH_M - VEH_F,
  sex_in_PAS = PAS_M - PAS_F,
  levels = design
)

PRIMARY <- c("PAS_vs_VEH_female", "PAS_vs_VEH_male", "treatment_by_sex")
SUPPLEMENTARY <- c("sex_in_VEH", "sex_in_PAS")

stopifnot(
  setequal(c(PRIMARY, SUPPLEMENTARY), colnames(contrasts)),
  all.equal(
    contrasts[, "treatment_by_sex"], contrasts[, "sex_in_VEH"] - contrasts[, "sex_in_PAS"]
  )
)

stopifnot(identical(rownames(design), colnames(proteins$E)))

universes <- list(full = rep(TRUE, nrow(proteins)), mito = proteins$genes$is_mito)

group_means <- t(apply(proteins$E, 1, \(x) tapply(x, proteins$targets$group, mean)))

manual <- cbind(
  PAS_vs_VEH_female = group_means[, "PAS_F"] - group_means[, "VEH_F"],
  PAS_vs_VEH_male = group_means[, "PAS_M"] - group_means[, "VEH_M"],
  treatment_by_sex = (group_means[, "PAS_F"] - group_means[, "VEH_F"]) -
    (group_means[, "PAS_M"] - group_means[, "VEH_M"]),
  sex_in_VEH = group_means[, "VEH_M"] - group_means[, "VEH_F"],
  sex_in_PAS = group_means[, "PAS_M"] - group_means[, "PAS_F"]
)

fitted <- contrasts.fit(lmFit(proteins$E, design), contrasts)$coefficients

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

saveRDS(
  list(
    design = design, contrasts = contrasts, universes = universes,
    primary = PRIMARY, supplementary = SUPPLEMENTARY
  ),
  file.path(DESIGN_DATA, "design.rds")
)

write_workbook(
  list(
    design = as_tibble(design, rownames = "sample_id"),
    contrasts = as_tibble(contrasts, rownames = "group"),
    contrast_check = contrast_check,
    roles = tibble(
      contrast = colnames(contrasts),
      role = if_else(colnames(contrasts) %in% PRIMARY, "primary", "supplementary")
    ),
    groups = count(proteins$targets, group, name = "n")
  ),
  c(
    "Cell-means design matrix, one row per mouse.",
    "Contrast matrix: the five questions as group weights.",
    "Each contrast checked against its intended group difference.",
    "Which contrasts are primary and which supplementary.",
    "Mice per group."
  ),
  file.path(DESIGN_DATA, "01_design.xlsx"),
  c(proteins = file.path(QUANT_DATA, "proteins.rds"))
)

contrast_check

sessionInfo()
