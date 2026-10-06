#!/usr/bin/env Rscript
# PAS against vehicle in males, phenotype links. R/arc.R holds the step; the README explains it.

source(here::here("R", "arc.R"))
arc_associations("male")

sessionInfo()
