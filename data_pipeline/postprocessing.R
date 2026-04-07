#!/usr/bin/env Rscript

# Reusable postprocessing helpers for downstream analyses

suppressPackageStartupMessages({
  library(tidyverse)
})

# Exclude lakes by name
# - df: data.frame with a LakeID column
# - lakes: character vector of LakeID values to exclude (case-insensitive)
exclude_lakes <- function(df, lakes = c("nimpkish", "msk"), .warn = TRUE) {
  if (!"LakeID" %in% names(df)) return(df)

  # Name-based exclusion
  extra_env <- Sys.getenv("EXCLUDE_LAKES_EXTRA", "")
  extra_opt <- getOption("exclude_lakes_extra", NULL)
  extra_vec <- character()
  if (nzchar(extra_env)) extra_vec <- c(extra_vec, strsplit(extra_env, ",")[[1]])
  if (!is.null(extra_opt)) extra_vec <- c(extra_vec, extra_opt)
  extra_vec <- tolower(trimws(extra_vec))
  lakes_all <- unique(c(tolower(lakes), extra_vec))
  remove_by_name <- tolower(df$LakeID) %in% lakes_all

  keep <- !remove_by_name
  out <- df[keep, , drop = FALSE]

  if (isTRUE(.warn)) {
    removed <- sum(!keep, na.rm = TRUE)
    removed_names <- unique(tolower(df$LakeID[remove_by_name]))
    msg <- paste0(
      "exclude_lakes: removed ", removed, " rows",
      if (length(removed_names)) paste0("; names: ", paste(sort(removed_names), collapse = ", ")) else ""
    )
    message(msg)
  }

  out
}

# Example usage downstream (not executed here):
# library(readr)
# source("data_pipeline/postprocessing.R")
# df <- readRDS("01_WorkingData/main_zooplankton_data.rds")
# df_clean <- exclude_lakes(df)                 # default excludes nimpkish and msk
# df_clean2 <- exclude_lakes(df, lakes = c("nimpkish","another"))
# # proceed with modeling using df_clean*


