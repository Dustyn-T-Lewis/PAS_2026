#!/usr/bin/env Rscript
# PAS against vehicle in females, protein by protein. R/arc.R holds the step; the README explains it.

source(here::here("R", "arc.R"))
arc_proteins("female")

sessionInfo()
