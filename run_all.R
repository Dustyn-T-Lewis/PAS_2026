# Runs the pipeline end to end, each script in its own R session so no step inherits objects or
# masked functions from the one before. Every script also runs alone; this file is a convenience.
# Folder numbers are not run order: the order below is.
#
#   Rscript run_all.R

steps <- c(
  "01_Preprocess/01_Filtering/a_script/01_filter.R",
  "01_Preprocess/02_Quantification/a_script/02_quantify.R",
  "02_Differential_Expression/01_Design/a_script/01_design.R",
  "02_Differential_Expression/02_Contrasts/a_script/02_contrasts.R",
  "03_Pathway_Enrichment/00_Gene_Sets/a_script/00_gene_sets.R",
  "03_Pathway_Enrichment/02_Contrasts/a_script/02_contrasts.R",
  "03_Pathway_Enrichment/01_Scores/a_script/01_scores.R",
  "03_Pathway_Enrichment/05_Volcano/a_script/05_volcano.R",
  "03_Pathway_Enrichment/06_NES_Scatter/a_script/06_nes_scatter.R",
  "03_Pathway_Enrichment/07_Compartments/a_script/07_compartments.R",
  "03_Pathway_Enrichment/08_Display/a_script/08_display.R",
  "04_Network/01_Build_Modules/a_script/01_build_modules.R",
  "04_Network/02_Contrasts/a_script/02_contrasts.R",
  "04_Network/05_Characterise/a_script/05_characterise.R",
  "04_Network/06_Mechanism/a_script/06_mechanism.R",
  "02_Differential_Expression/04_Associate/a_script/04_associate.R",
  "03_Pathway_Enrichment/04_Associate/a_script/04_associate.R",
  "04_Network/04_Associate/a_script/04_associate.R",
  "02_Differential_Expression/03_Classify/a_script/03_classify.R",
  "03_Pathway_Enrichment/03_Classify/a_script/03_classify.R",
  "04_Network/03_Classify/a_script/03_classify.R",
  "05_Summary/01_Signal/a_script/01_signal.R",
  "05_Summary/02_Panels/a_script/02_panels.R",
  "05_Summary/03_Biclusters/a_script/03_biclusters.R",
  "05_Summary/04_Sex_Patterns/a_script/04_sex_patterns.R",
  "05_Summary/05_Global/a_script/05_global.R",
  "05_Summary/06_Magnitude/a_script/06_magnitude.R",
  "05_Summary/07_Signature/a_script/07_signature.R",
  "05_Summary/08_Signature_Phenotype/a_script/08_signature_phenotype.R",
  file.path("06_Figures/F01/a_script", c("main/F01.R", "supp/S1.R", "supp/S2.R", "supp/S3.R")),
  file.path("06_Figures/F02/a_script", c("main/F02.R", "supp/S4.R")),
  "06_Figures/F03/a_script/main/F03.R",
  "06_Figures/F04/a_script/main/F04.R",
  "06_Figures/F05/a_script/main/F05.R"
)

root <- here::here()
for (step in steps) {
  message(format(Sys.time(), "%H:%M:%S "), step)
  status <- system2("Rscript", file.path(root, step), stdout = FALSE, stderr = FALSE)
  if (status != 0) stop("failed: ", step, " (rerun it alone to see the error)")
}
message("done")
