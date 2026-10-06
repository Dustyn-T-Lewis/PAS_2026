#!/usr/bin/env Rscript
# What PAS did to the phenotypes in females, and which proteins, pathways and protein groups move
# with them. R/arc.R holds the step; the README explains it.

source(here::here("R", "arc.R"))
arc_associations("female")

sessionInfo()
