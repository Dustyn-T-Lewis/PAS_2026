# Whether the PAS response differs by sex: (PAS_F − VEH_F) − (PAS_M − VEH_M).
# Runs the steps written once in 03_Comparisons/code/ for this comparison; phenotypes and
# classification do not apply here (see the README).
#   Rscript 03_Comparisons/Sex_Difference/a_script/run_comparison.R

source(here::here("03_Comparisons", "code", "01_proteins.R"))
source(here::here("03_Comparisons", "code", "02_pathways.R"))

report_proteins("sex_difference")
report_pathways("sex_difference")
