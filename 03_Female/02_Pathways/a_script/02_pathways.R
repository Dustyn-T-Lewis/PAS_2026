#!/usr/bin/env Rscript
# Which pathways PAS shifted in females.
# The step is report_pathways() in R/comparisons.R.

source(here::here("R", "comparisons.R"))
report_pathways("female")

sessionInfo()
