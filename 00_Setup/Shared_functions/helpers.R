# Definitions and writers every stage shares, so a contrast, a phenotype or a workbook layout is
# defined once. Every script sources this file; it sources figures.R.

library(here)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(tidyr)
source(here::here("00_Setup", "Shared_functions", "figures.R"))

SEED <- 20260910
ALPHA <- 0.05
CONTRASTS <- c(
  "PAS_vs_VEH_female", "PAS_vs_VEH_male", "treatment_by_sex", "sex_in_VEH", "sex_in_PAS"
)
# Each comparison reads one contrast of the shared fit; the two baseline sex contrasts are context.
COMPARISONS <- list(
  female = list(dir = "Female", contrast = "PAS_vs_VEH_female", sex = "F", label = "Female"),
  male = list(dir = "Male", contrast = "PAS_vs_VEH_male", sex = "M", label = "Male"),
  sex_difference = list(
    dir = "Sex_Difference", contrast = "treatment_by_sex", sex = NA, label = "Sex difference"
  )
)
COMPARISON_CONTRASTS <- map_chr(COMPARISONS, "contrast")
COMPARISON_OF_CONTRAST <- c(
  PAS_vs_VEH_female = "female", PAS_vs_VEH_male = "male", treatment_by_sex = "sex_difference",
  sex_in_VEH = "context", sex_in_PAS = "context"
)
PHENOTYPES <- c(
  "capillary_density", "csa_mean", "satellite_cells_per_fiber",
  "conductance_sol_pyr", "conductance_pla_pyr", "conductance_sol_oct"
)

input_path <- \(...) here("00_Setup", "Input", ...)
output_path <- \(stage, step, file) here(stage, step, "c_data", file)

stage_paths <- function(stage, step, reports = TRUE) {
  root <- here(stage, step)
  dirs <- file.path(root, c("c_data", if (reports) "b_reports"))
  walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  list(
    root = root, data = file.path(root, "c_data"), reports = file.path(root, "b_reports"),
    workbook = file.path(root, "c_data", paste0(tolower(step), ".xlsx"))
  )
}

# A set is tested when this many of its genes are detected in the universe (enrichVolcano's
# min_size and max_size); 04_Gene_Sets trims the lists to it and 02_Enrichment passes it to fgsea.
TEST_WINDOW <- list(full = c(15, 500), mito = c(10, 500))

# DIA-NN lists every gene of a protein group; the first names the group.
first_symbol <- \(genes) str_split_i(genes, ";", 1)

# Set names in MSigDB style: upper case, runs of other characters as one underscore.
set_id <- \(x) str_to_upper(str_replace_all(x, "[^A-Za-z0-9]+", "_"))

# STRING v12, mouse, medium confidence, read from the copies in 00_Setup/Input/STRING_v12.
STRING_FILES <- set_names(
  list.files(input_path("STRING_v12"), full.names = TRUE),
  c("string_aliases", "string_info", "string_links")
)
string_db <- function() {
  STRINGdb::STRINGdb$new(
    version = "12.0", species = 10090, score_threshold = 400, network_type = "full",
    input_directory = input_path("STRING_v12")
  )
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
      str_remove(str_to_upper(label), "^(HALLMARK|REACTOME)_"),
      40,
      ellipsis = ""
    ), .after = label)
}

# Every workbook opens on read_me and ends with the inputs' md5 and the packages attached.
write_workbook <- function(sheets, holds, file, inputs) {
  stopifnot(length(sheets) == length(holds))
  sheets <- c(sheets, list(
    input_manifest = tibble(
      input = names(inputs), path = str_remove(unname(inputs), fixed(paste0(here(), "/"))),
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
