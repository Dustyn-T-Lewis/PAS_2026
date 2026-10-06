# Definitions and writers every stage shares, so a contrast, a phenotype or a workbook layout is
# defined once. Sourced as source(here::here("R", "helpers.R")); it sources R/figures.R.

library(here)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(tidyr)
source(here::here("R", "figures.R"))

SEED <- 20260910
ALPHA <- 0.05
GROUPS <- c("VEH_F", "VEH_M", "PAS_F", "PAS_M")
CONTRASTS <- c(
  "PAS_vs_VEH_female", "PAS_vs_VEH_male", "treatment_by_sex", "sex_in_VEH", "sex_in_PAS"
)
# Each arc reads one contrast of the shared fit; the two baseline sex contrasts are trunk context.
ARCS <- list(
  female = list(dir = "01_Female", contrast = "PAS_vs_VEH_female", sex = "F", label = "Female"),
  male = list(dir = "02_Male", contrast = "PAS_vs_VEH_male", sex = "M", label = "Male"),
  sex_difference = list(
    dir = "03_Sex_Difference", contrast = "treatment_by_sex", sex = NA, label = "Sex difference"
  )
)
ARC_OF_CONTRAST <- c(
  PAS_vs_VEH_female = "female", PAS_vs_VEH_male = "male", treatment_by_sex = "sex_difference",
  sex_in_VEH = "trunk", sex_in_PAS = "trunk"
)
# Two-group comparisons for r with group: first group = 1.
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

input_path <- \(...) here("00_Trunk", "00_Input", ...)

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
  ph <- readxl::read_excel(input_path("phenotypes.xlsx"), "phenotypes")
  ph <- ph[match(samples, ph$mouse), ]
  stopifnot(identical(ph$mouse, samples), !anyNA(ph[PHENOTYPES]))
  ph
}

# Π (Xiao et al. 2014): p raised to |log2 fold change|. A ranking that rewards large changes;
# it controls no error rate, and Π <= 0.05 can include proteins with p > 0.05.
pi_score <- \(p, logfc) p^abs(logfc)

tier_counts <- function(p, fdr, pi = NULL) {
  tibble(
    tests = length(p), p_le_0.05 = sum(p <= ALPHA), expected_by_chance = ALPHA * length(p),
    pi_le_0.05 = if (is.null(pi)) NA_integer_ else sum(pi <= ALPHA),
    fdr_le_0.10 = sum(fdr <= 0.10), fdr_le_0.05 = sum(fdr <= ALPHA)
  )
}

# One protein per gene, the rule enrichVolcano applies: the most abundant protein represents the
# gene, ties broken by accession. Returns a logical keep vector in input order.
one_per_gene <- function(protein, gene, abundance) {
  o <- order(gene, -abundance, protein)
  keep <- logical(length(protein))
  keep[o] <- !duplicated(gene[o])
  keep & !is.na(gene) & gene != ""
}

# Reactome's disease and infection branch, left out of the tested Reactome sets and of labels.
REACTOME_DISEASE <- paste0(
  "INFECTION|INFECTIOUS|DISEASE|VIRUS|VIRAL|CANCER|BACTERIAL|LEISHMANIA|SARS_COV|HIV|",
  "INFLUENZA|DEFECTIVE|MUTANT|VARIANT"
)

# A protein group's name: its most over-represented term (fgsea::fora against the detected
# proteome, BH across the sets tested for that group) when that term reaches FDR <= 0.05, else
# "unannotated". Used for STRING clusters and modules, so no group is named by a chance hit.
# Reactome's disease and infection branch is left out, as in the figures.
label_groups <- function(groups, collections, universe) {
  sets <- unlist(unname(collections), recursive = FALSE)
  sets <- sets[!str_detect(names(sets), paste0("^REACTOME_.*(", REACTOME_DISEASE, ")"))]
  map(groups, \(g) {
    top <- fgsea::fora(sets, g, universe, minSize = 10, maxSize = 500) |>
      as_tibble() |>
      slice_min(pval, n = 1, with_ties = FALSE)
    named <- nrow(top) == 1 && top$padj <= ALPHA
    tibble(
      size = length(g), label = if (named) top$pathway else "unannotated",
      label_fdr = if (named) top$padj else NA_real_, members = paste(g, collapse = ";")
    )
  }) |>
    list_rbind() |>
    mutate(label_id = str_trunc(
      str_remove(str_to_upper(label), "^(HALLMARK|REACTOME|GOBP|KEGG_MEDICUS_REFERENCE)_"),
      40,
      ellipsis = ""
    ), .after = label)
}

# Rows ordered by hierarchical clustering on 1 - Pearson r, average linkage: display only.
cluster_order <- \(m) hclust(as.dist(1 - cor(t(m), use = "pairwise.complete.obs")), "average")$order

# Every workbook opens on read_me and ends with the inputs' md5 and the packages loaded.
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
