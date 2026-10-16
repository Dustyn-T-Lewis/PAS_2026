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
# DIA-NN precursor and library q-value cutoff, applied in 01_Filtering and again when
# 05_Composition reads the report for haemoglobin.
Q_CUTOFF <- 0.01
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
SEX_COMPARISONS <- COMPARISONS[c("female", "male")]
COMPARISON_OF_CONTRAST <- c(
  set_names(names(COMPARISONS), COMPARISON_CONTRASTS),
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

# dpcQuant takes about two and a half minutes per detection slope, so each fit is kept in c_data/
# and refitted when the filtered precursors change or REFIT_DPCQUANT=TRUE.
dpcquant_cache <- \(dir, slope) file.path(dir, sprintf("dpcquant_slope_%s.rds", slope))
dpcquant_cached <- function(y, slope, file) {
  fingerprint <- c(dim(y), sum(y$E, na.rm = TRUE))
  refit <- isTRUE(as.logical(Sys.getenv("REFIT_DPCQUANT", "FALSE")))
  cached <- if (file.exists(file) && !refit) readRDS(file)
  if (is.null(cached) || !identical(attr(cached, "fingerprint"), fingerprint)) {
    cached <- limpa::dpcQuant(y, "Protein.Group", dpc.slope = slope, verbose = FALSE)
    attr(cached, "fingerprint") <- fingerprint
    saveRDS(cached, file)
  }
  cached
}

# Median centring is written out: limma's "scale" divides by column medians, which suits unlogged
# intensities, not log2 values.
normalize <- function(E, method) {
  switch(method,
    none = E,
    median = {
      medians <- apply(E, 2, median)
      sweep(E, 2, medians) + mean(medians)
    },
    limma::normalizeBetweenArrays(E, method)
  )
}

# Run depth, the covariate the surrogate variable tracks: log2 precursors identified in each run.
run_depth_log2 <- function(run_qc, sample_id) {
  log2(run_qc$precursors[match(sample_id, run_qc$sample_id)])
}

# DIA-NN lists every gene of a protein group; the first names the group.
first_symbol <- \(genes) str_split_i(genes, ";", 1)

# Set names in MSigDB style: upper case, runs of other characters as one underscore.
set_id <- \(x) str_to_upper(str_replace_all(x, "[^A-Za-z0-9]+", "_"))

# STRING v12, mouse, medium confidence, read from the copies in 00_Setup/Input/STRING_v12.
STRING_FILES <- map_chr(c(
  string_aliases = "10090.protein.aliases.v12.0.txt.gz",
  string_info = "10090.protein.info.v12.0.txt.gz",
  string_links = "10090.protein.links.v12.0.min400.txt.gz"
), \(f) input_path("STRING_v12", f))
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
LABEL_MIN_SIZE <- 10
label_groups <- function(groups, collections, universe) {
  sets <- unlist(unname(collections), recursive = FALSE)
  sets <- sets[!str_detect(names(sets), paste0("^REACTOME_.*(", REACTOME_DISEASE, ")"))]
  map(groups, \(g) {
    top <- fgsea::fora(sets, g, universe, minSize = LABEL_MIN_SIZE, maxSize = 500) |>
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

# A step's preview PDF: one landscape page per named function, each drawing with a package's own
# plot or base graphics; cairo, so Π and ≤ render. Returns the rows of the workbook's plot_index.
preview_pdf <- function(file, pages) {
  cairo_pdf(file, width = 7, height = 5, onefile = TRUE)
  on.exit(dev.off())
  defaults <- par(no.readonly = TRUE)
  walk(pages, \(draw) {
    draw()
    par(defaults)
  })
  tibble(file = basename(file), page = seq_along(pages), shows = names(pages))
}

# fgsea's own table of running scores for the n strongest sets each way among results (term, NES,
# p, padj, direction), on the ranking they were tested against; a preview_pdf page.
gsea_table_page <- function(results, sets, ranks, n = 5) {
  top <- results |>
    slice_min(p, n = n, by = direction, with_ties = FALSE) |>
    arrange(desc(NES))
  short <- set_label(top$term, 36)
  print(fgsea::plotGseaTable(
    set_names(sets[top$term], short), ranks,
    transmute(top, pathway = short, NES, pval = p, padj),
    gseaParam = 1
  ))
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
