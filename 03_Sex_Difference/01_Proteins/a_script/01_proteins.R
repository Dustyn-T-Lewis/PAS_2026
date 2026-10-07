#!/usr/bin/env Rscript
# Which proteins respond to PAS differently in females and males.
# The step is arc_proteins() in R/arc.R.

source(here::here("R", "arc.R"))
arc_proteins("sex_difference")

sessionInfo()
