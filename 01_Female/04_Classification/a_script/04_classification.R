#!/usr/bin/env Rscript
# Which single features separate PAS from vehicle in females. R/arc.R holds the step; the README
# explains it.

source(here::here("R", "arc.R"))
arc_classification("female")

sessionInfo()
