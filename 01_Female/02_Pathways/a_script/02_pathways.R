#!/usr/bin/env Rscript
# Which pathways PAS shifted in females.
# The step is arc_pathways() in R/arc.R.

source(here::here("R", "arc.R"))
arc_pathways("female")

sessionInfo()
