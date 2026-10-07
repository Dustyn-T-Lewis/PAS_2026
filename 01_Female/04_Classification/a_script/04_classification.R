#!/usr/bin/env Rscript
# Which single scores separate PAS from vehicle in females.
# The step is arc_classification() in R/arc.R.

source(here::here("R", "arc.R"))
arc_classification("female")

sessionInfo()
