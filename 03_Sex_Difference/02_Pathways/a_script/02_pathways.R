#!/usr/bin/env Rscript
# Sex difference in the PAS response (female minus male), pathway by pathway. R/arc.R holds the
# step; the README explains it.

source(here::here("R", "arc.R"))
arc_pathways("sex_difference")

sessionInfo()
