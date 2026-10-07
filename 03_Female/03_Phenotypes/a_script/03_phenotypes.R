#!/usr/bin/env Rscript
# Phenotype effects of PAS in females, and the proteins and scores associated with each phenotype.
# The step is report_phenotypes() in R/comparisons.R.

source(here::here("R", "comparisons.R"))
report_phenotypes("female")

sessionInfo()
