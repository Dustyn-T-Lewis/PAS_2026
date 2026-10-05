# Definitions and writers every stage shares, so a contrast, a phenotype or a workbook layout is
# defined once. Sourced as source(here::here("R", "helpers.R")).

library(here)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(tidyr)

SEED <- 20260910
ALPHA <- 0.05
GROUPS <- c("VEH_F", "VEH_M", "PAS_F", "PAS_M")
CONTRASTS <- c(
  "PAS_vs_VEH_female", "PAS_vs_VEH_male", "treatment_by_sex", "sex_in_VEH", "sex_in_PAS"
)
PRIMARY <- CONTRASTS[1:3]
# Two-group comparisons for classification and correlation tests: first group minus second.
COMPARISONS <- list(
  PAS_vs_VEH_female = c("PAS_F", "VEH_F"),
  PAS_vs_VEH_male = c("PAS_M", "VEH_M"),
  sex_in_VEH = c("VEH_M", "VEH_F"),
  sex_in_PAS = c("PAS_M", "PAS_F")
)
PHENOTYPES <- c(
  "capillary_density", "csa_mean", "satellite_cells_per_fiber",
  "conductance_sol_pyr", "conductance_pla_pyr", "conductance_sol_oct"
)

stage_paths <- function(stage, step) {
  root <- here(stage, step)
  dirs <- file.path(root, c("b_reports", "c_data"))
  walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  list(
    root = root, data = file.path(root, "c_data"), reports = file.path(root, "b_reports"),
    workbook = file.path(root, "c_data", paste0(tolower(step), ".xlsx"))
  )
}

read_phenotypes <- function(samples) {
  ph <- readxl::read_excel(here("00_Input", "phenotypes.xlsx"), "phenotypes")
  ph <- ph[match(samples, ph$mouse), ]
  stopifnot(identical(ph$mouse, samples), !anyNA(ph[PHENOTYPES]))
  ph
}

# Π (Xiao et al. 2014): p raised to |log2 fold change|. A ranking that rewards large changes;
# it controls no error rate.
pi_score <- \(p, logfc) p^abs(logfc)

tier_counts <- function(p, fdr, pi = NULL) {
  tibble(
    tests = length(p), p_le_0.05 = sum(p <= ALPHA), expected_by_chance = ALPHA * length(p),
    pi_le_0.05 = if (is.null(pi)) NA_integer_ else sum(pi <= ALPHA),
    fdr_le_0.10 = sum(fdr <= 0.10), fdr_le_0.05 = sum(fdr <= ALPHA)
  )
}

# Labels from each sex's effect and the direct interaction test. Only the interaction claims the
# sexes differ; a hit in one sex alone is "detected in" that sex.
sex_label <- function(effect_f, p_f, effect_m, p_m, p_interaction) {
  case_when(
    p_interaction <= ALPHA ~ "sex-differential",
    p_f <= ALPHA & p_m <= ALPHA & sign(effect_f) == sign(effect_m) ~ "shared",
    p_f <= ALPHA & p_m > ALPHA ~ "detected in females",
    p_m <= ALPHA & p_f > ALPHA ~ "detected in males",
    p_f <= ALPHA & p_m <= ALPHA ~ "opposite, interaction not significant",
    TRUE ~ "neither"
  )
}

# Phenotype designs on top of the cell-means design: one shared slope; the slope's difference in
# males; or the female mice alone.
with_phenotype <- \(design, x) cbind(design, phenotype = as.numeric(scale(x)))
with_sex_slope <- function(design, x, male) {
  z <- as.numeric(scale(x))
  cbind(design, phenotype = z, phenotype_male_minus_female = z * male)
}

# One protein per gene, the rule enrichVolcano applies: the most abundant protein represents the
# gene, ties broken by accession. Returns a logical keep vector in input order.
one_per_gene <- function(protein, gene, abundance) {
  o <- order(gene, -abundance, protein)
  keep <- logical(length(protein))
  keep[o] <- !duplicated(gene[o])
  keep & !is.na(gene) & gene != ""
}

muffle <- function(expr, patterns) {
  withCallingHandlers(expr, warning = \(w) {
    if (any(str_detect(conditionMessage(w), patterns))) invokeRestart("muffleWarning")
  })
}

# Rows ordered by hierarchical clustering on 1 - Pearson r, average linkage: display only.
cluster_order <- \(m) hclust(as.dist(1 - cor(t(m), use = "pairwise.complete.obs")), "average")$order

overview_heatmap <- function(m, group, title) {
  z <- t(scale(t(m)))
  z <- z[cluster_order(z), order(factor(group, GROUPS)), drop = FALSE]
  limit <- unname(quantile(abs(z), 0.99, na.rm = TRUE))
  ComplexHeatmap::Heatmap(z,
    name = "z", col = circlize::colorRamp2(c(-limit, 0, limit), c("#2166AC", "white", "#B2182B")),
    cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = nrow(z) <= 60,
    column_split = factor(group[order(factor(group, GROUPS))], GROUPS), column_title_gp =
      grid::gpar(fontsize = 8), column_names_gp = grid::gpar(fontsize = 6),
    row_names_gp = grid::gpar(fontsize = 5), use_raster = nrow(z) > 200, column_title = title
  )
}

# One PDF from a named list of drawing functions, one page each; returns the rows plot_index needs.
# cairo_pdf, because the base device cannot draw Π.
write_pdf <- function(file, pages, width = 11, height = 8.5) {
  cairo_pdf(file, width = width, height = height, onefile = TRUE)
  on.exit(dev.off())
  iwalk(pages, \(draw, title) draw())
  tibble(file = basename(file), page = seq_along(pages), shows = names(pages))
}

# Every workbook opens on read_me and ends with the inputs' md5 and the packages loaded.
# Terms drawn in summary figures: Hallmark, Reactome without its disease and infection branch,
# KEGG MEDICUS reference pathways, GO:BP without the broad GO Slim terms, and every MitoCarta
# pathway. Display only; tests and corrections cover every term.
REACTOME_DISEASE <- paste0(
  "INFECTION|INFECTIOUS|DISEASE|VIRUS|VIRAL|CANCER|BACTERIAL|LEISHMANIA|SARS_COV|HIV|",
  "INFLUENZA|DEFECTIVE|MUTANT|VARIANT"
)
display_pool <- function(results, slim_terms) {
  filter(
    results,
    database %in% c("Hallmark", "MitoCarta") |
      (database == "Reactome" & !str_detect(term, REACTOME_DISEASE)) |
      (database == "KEGG" & str_starts(term, "KEGG_MEDICUS_REFERENCE")) |
      (database == "GO:BP" & !term %in% slim_terms)
  )
}

write_workbook <- function(sheets, holds, file, inputs) {
  stopifnot(length(sheets) == length(holds))
  sheets <- c(sheets, list(
    input_manifest = tibble(
      input = names(inputs), path = str_remove(unname(inputs), paste0(here(), "/")),
      md5 = unname(tools::md5sum(inputs))
    ),
    package_versions = sessioninfo::package_info("attached", dependencies = FALSE) |>
      as_tibble() |>
      dplyr::select(package, version = loadedversion, source)
  ))
  read_me <- tibble(
    sheet = names(sheets),
    holds = c(holds, "Files read, with md5.", "Packages attached at run time.")
  )
  writexl::write_xlsx(c(list(read_me = read_me), sheets), file)
}

# Each feature's Pearson r with one phenotype inside each group of a comparison, the two compared
# by Fisher's z (DGCA). The permutation q-value needs many features, so small sets pass "BH".
differential_correlation <- function(features, ph, adjust = "perm", n_perm = 200) {
  expand_grid(comparison = names(COMPARISONS), phenotype = PHENOTYPES) |>
    pmap(\(comparison, phenotype) {
      groups <- COMPARISONS[[comparison]]
      keep <- ph$group %in% groups
      design <- model.matrix(~ 0 + factor(ph$group[keep], levels = groups))
      colnames(design) <- groups
      res <- suppressMessages(DGCA::ddcorAll(
        inputMat = rbind(features[, keep, drop = FALSE], phenotype = ph[[phenotype]][keep]),
        design = design, compare = groups, splitSet = "phenotype", corrType = "pearson",
        adjust = adjust, nPerm = if (adjust == "perm") n_perm else 0, verbose = FALSE
      ))
      tibble(
        comparison = comparison, phenotype = phenotype, feature = res$Gene1,
        r_first = res[[paste0(groups[1], "_cor")]], r_second = res[[paste0(groups[2], "_cor")]],
        z_diff = res$zScoreDiff, p_diff = res$pValDiff, q_diff = res$pValDiff_adj,
        class = res$Classes
      )
    }) |>
    list_rbind()
}

# How well each feature alone separates the two groups of each comparison: AUC with DeLong's
# interval (pROC; undefined at perfect separation) and the exact Wilcoxon p, BH within comparison.
# At five against five, chance alone separates a feature perfectly with probability 2/252.
auc_screen <- function(m, group) {
  imap(COMPARISONS, \(groups, cmp) {
    keep <- group %in% groups
    y <- factor(group[keep], levels = rev(groups))
    map(rownames(m), \(f) {
      x <- m[f, keep]
      roc <- pROC::roc(y, x, levels = levels(y), direction = "<", quiet = TRUE)
      ci <- muffle(as.numeric(pROC::ci.auc(roc, method = "delong")), "ci.auc|DeLong|NaN")
      tibble(
        feature = f, auc = as.numeric(roc$auc), auc_low = ci[1], auc_high = ci[3],
        p = wilcox.test(x[y == groups[1]], x[y == groups[2]], exact = TRUE)$p.value
      )
    }) |>
      list_rbind() |>
      mutate(comparison = cmp, fdr = p.adjust(p, "BH"), .before = 1)
  }) |>
    list_rbind()
}

auc_summary <- function(auc) {
  auc |>
    summarise(
      features = n(), p_le_0.05 = sum(p <= ALPHA), fdr_le_0.10 = sum(fdr <= 0.10),
      perfect = sum(auc %in% c(0, 1)), perfect_expected_by_chance = n() * 2 / choose(10, 5),
      .by = comparison
    )
}

# Strip charts of the strongest separators, one page per comparison.
auc_pages <- function(auc, m, group, label = identity, n = 12) {
  imap(COMPARISONS, \(groups, cmp) {
    top <- auc |>
      filter(comparison == cmp) |>
      slice_min(p, n = n, with_ties = FALSE)
    \() {
      par(mfrow = c(3, 4), mar = c(3, 4, 3, 1))
      pwalk(top, \(feature, auc, p, ...) {
        keep <- group %in% groups
        stripchart(m[feature, keep] ~ factor(group[keep], groups),
          vertical = TRUE, method = "jitter",
          pch = 19, cex = 0.8, ylab = "",
          main = sprintf("%s\nAUC %.2f, p %.3f", str_trunc(label(feature), 30), auc, p),
          cex.main = 0.7
        )
      })
    }
  }) |>
    set_names(str_glue("12 strongest separators, {names(COMPARISONS)}"))
}
