# Runs the pipeline end to end, each script in its own R session so no step inherits objects or
# masked functions from the one before. Every script also runs alone; this file is a convenience.
# Processing and statistics run first; the three comparisons report their results; the figures
# read the statistics. Output from every script goes to run_all.log.
#
#   Rscript run_all.R

processing <- file.path("01_Processing", c(
  "01_Filtering/a_script/01_filtering.R",
  "02_Quantification/a_script/02_quantification.R",
  "03_Design/a_script/03_design.R",
  "04_Gene_Sets/a_script/04_gene_sets.R",
  "05_Scores/a_script/05_scores.R",
  "06_Modules/a_script/06_modules.R"
))
statistics <- file.path("02_Statistics", c(
  "01_Fit/a_script/01_fit.R",
  "02_Enrichment/a_script/02_enrichment.R",
  "03_Classification/a_script/03_classification.R",
  "04_PCA/a_script/04_pca.R"
))
comparisons <- file.path(
  "03_Comparisons", c("Female", "Male", "Sex_Difference"), "a_script", "run_comparison.R"
)
figures <- "04_Figures/F01_Overview/a_script/f01_overview.R"

root <- here::here()
log_file <- file.path(root, "run_all.log")
invisible(file.create(log_file))
for (step in c(processing, statistics, comparisons, figures)) {
  message(format(Sys.time(), "%H:%M:%S "), step)
  status <- system2("Rscript", c(shQuote(file.path(root, step)), ">>", shQuote(log_file), "2>&1"))
  if (status != 0) stop("failed: ", step, " (see run_all.log)")
}
message("done")
