#!/usr/bin/env Rscript
# Phenotype effects of PAS in females, and the proteins and scores associated with each phenotype.
# The step is arc_associations() in R/arc.R.

source(here::here("R", "arc.R"))
arc_associations("female")

sessionInfo()
