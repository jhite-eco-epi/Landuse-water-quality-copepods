#!/usr/bin/env Rscript

source("correlations/helpers_correlations.R")

if (sys.nframe() == 0) {
  message("Building all correlation outputs...")
  result <- build_all_correlations()
  message("Results README written to: ", result$readme)
}
