#!/usr/bin/env Rscript

source("correlations/helpers_correlations.R")

if (sys.nframe() == 0) {
  message("Building abiotic correlation outputs...")
  build_abiotic_correlations()
}
