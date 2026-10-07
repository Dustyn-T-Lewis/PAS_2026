#!/usr/bin/env Rscript
# Which proteins PAS changed in females.
# The step is report_proteins() in R/comparisons.R.

source(here::here("R", "comparisons.R"))
report_proteins("female")

sessionInfo()
