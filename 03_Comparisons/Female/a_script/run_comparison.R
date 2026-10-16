# PAS against vehicle in females (PAS_F − VEH_F).
# Runs the steps written once in 03_Comparisons/code/ for this comparison.
#   Rscript 03_Comparisons/Female/a_script/run_comparison.R

source(here::here("03_Comparisons", "code", "01_proteins.R"))
source(here::here("03_Comparisons", "code", "02_pathways.R"))
source(here::here("03_Comparisons", "code", "03_phenotypes.R"))
source(here::here("03_Comparisons", "code", "04_classification.R"))

report_proteins("female")
report_pathways("female")
report_phenotypes("female")
report_classification("female")
