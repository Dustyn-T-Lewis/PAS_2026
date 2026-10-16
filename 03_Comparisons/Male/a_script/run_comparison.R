# PAS against vehicle in males (PAS_M − VEH_M).
# Runs the steps written once in 03_Comparisons/code/ for this comparison.
#   Rscript 03_Comparisons/Male/a_script/run_comparison.R

source(here::here("03_Comparisons", "code", "01_proteins.R"))
source(here::here("03_Comparisons", "code", "02_pathways.R"))
source(here::here("03_Comparisons", "code", "03_phenotypes.R"))
source(here::here("03_Comparisons", "code", "04_classification.R"))

report_proteins("male")
report_pathways("male")
report_phenotypes("male")
report_classification("male")
