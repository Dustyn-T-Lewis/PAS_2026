# Runs the pipeline end to end, each script in its own R session so no step inherits objects or
# masked functions from the one before. Every script also runs alone; this file is a convenience.
# The trunk runs first and computes everything; the arcs slice its results; the figures read
# both.
#
#   Rscript run_all.R

trunk <- file.path("00_Trunk", c(
  "01_Filtering/a_script/01_filtering.R",
  "02_Quantification/a_script/02_quantification.R",
  "03_Design/a_script/03_design.R",
  "04_Gene_Sets/a_script/04_gene_sets.R",
  "05_Scores/a_script/05_scores.R",
  "06_Modules/a_script/06_modules.R",
  "07_Fit/a_script/07_fit.R",
  "08_Enrichment/a_script/08_enrichment.R",
  "09_Classification/a_script/09_classification.R",
  "10_Global/a_script/10_global.R"
))
arc_steps <- c(
  "01_Proteins/a_script/01_proteins.R",
  "02_Pathways/a_script/02_pathways.R",
  "03_Associations/a_script/03_associations.R",
  "04_Classification/a_script/04_classification.R"
)
arcs <- c(
  file.path("01_Female", arc_steps),
  file.path("02_Male", arc_steps),
  file.path("03_Sex_Difference", arc_steps[1:2])
)
figures <- "Figures/F01_Overview/a_script/f01_overview.R"

root <- here::here()
for (step in c(trunk, arcs, figures)) {
  message(format(Sys.time(), "%H:%M:%S "), step)
  status <- system2("Rscript", file.path(root, step), stdout = FALSE, stderr = FALSE)
  if (status != 0) stop("failed: ", step, " (rerun it alone to see the error)")
}
message("done")
