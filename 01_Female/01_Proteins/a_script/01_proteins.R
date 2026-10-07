#!/usr/bin/env Rscript
# Which proteins PAS changed in females.
# The step is arc_proteins() in R/arc.R.

source(here::here("R", "arc.R"))
arc_proteins("female")

sessionInfo()
