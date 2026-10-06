#!/usr/bin/env Rscript
# PAS against vehicle in males, single-feature classification. R/arc.R holds the step; the README explains it.

source(here::here("R", "arc.R"))
arc_classification("male")

sessionInfo()
