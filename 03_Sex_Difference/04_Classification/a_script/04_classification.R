#!/usr/bin/env Rscript
# Which single features separate males from females within vehicle, within PAS and on each PAS
# mouse's change. R/arc.R holds the step; the README explains it.

source(here::here("R", "arc.R"))
arc_classification("sex_difference")

sessionInfo()
