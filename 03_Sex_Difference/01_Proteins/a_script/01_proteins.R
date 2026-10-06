#!/usr/bin/env Rscript
# Sex difference in the PAS response (female minus male), protein by protein. R/arc.R holds the
# step; the README explains it.

source(here::here("R", "arc.R"))
arc_proteins("sex_difference")

sessionInfo()
