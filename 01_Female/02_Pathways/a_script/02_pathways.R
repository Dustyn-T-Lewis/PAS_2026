#!/usr/bin/env Rscript
# PAS against vehicle in females, pathway by pathway. R/arc.R holds the step; the README explains it.

source(here::here("R", "arc.R"))
arc_pathways("female")

sessionInfo()
