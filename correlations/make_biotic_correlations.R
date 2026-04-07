#!/usr/bin/env Rscript

source("correlations/helpers_correlations.R")

if (sys.nframe() == 0) {
  message("Building biotic correlation outputs...")
  build_biotic_correlations()
}
