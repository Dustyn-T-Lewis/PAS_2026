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
  "04_PCA/a_script/04_pca.R",
  "05_Composition/a_script/05_composition.R",
  "06_Candidates/a_script/06_candidates.R",
  "07_Sensitivity/a_script/07_sensitivity.R"
))
comparisons <- file.path(
  "03_Comparisons", c("Female", "Male", "Sex_Difference"), "a_script", "run_comparison.R"
)
figures <- file.path("04_Figures", c(
  "F01_Global/a_script/F01_composite.R",
  "F02_Proteins/a_script/F02_composite.R",
  "F03_Pathways/a_script/F03_composite.R",
  "F04_Phenotype_pathways/a_script/F04_composite.R",
  "F05_Modules/a_script/F05_composite.R",
  "F05_Complexes/a_script/F05_composite.R",
  "S1_Processing/a_script/S1_composite.R",
  "S2_Proteins/a_script/S2_composite.R",
  "S3_Pathway_sets/a_script/S3_composite.R",
  "S4_Protein_associations/a_script/S4_composite.R",
  "S5_Structure/a_script/S5_composite.R",
  "S6_Protein_patterns/a_script/S6_composite.R",
  "S7_Candidates/a_script/S7_composite.R",
  "S8_Sensitivity/a_script/S8_composite.R"
))

root <- here::here()
log_file <- file.path(root, "run_all.log")
invisible(file.create(log_file))
for (step in c(processing, statistics, comparisons, figures)) {
  message(format(Sys.time(), "%H:%M:%S "), step)
  status <- system2("Rscript", c(shQuote(file.path(root, step)), ">>", shQuote(log_file), "2>&1"))
  if (status != 0) stop("failed: ", step, " (see run_all.log)")
}
message("done")
