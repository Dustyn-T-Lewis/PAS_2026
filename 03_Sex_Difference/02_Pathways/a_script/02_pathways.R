#!/usr/bin/env Rscript
# Which pathways respond to PAS differently in females and males.
# The step is arc_pathways() in R/arc.R.

source(here::here("R", "arc.R"))
arc_pathways("sex_difference")

sessionInfo()
