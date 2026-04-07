#!/usr/bin/env Rscript

# Reusable postprocessing helpers for downstream analyses

suppressPackageStartupMessages({
  library(tidyverse)
})

# Shared default exclusions for downstream analyses and reporting
default_lake_exclusions <- function() {
  tibble::tibble(
    LakeID = c("Nimpkish", "Msk"),
    reason = c(
      "Saline/ultra-deep (a priori exclusion)",
      "Manual exclusion"
    )
  )
}

# Exclude lakes by name
# - df: data.frame with a lake identifier column
# - lakes: character vector of LakeID values to exclude (case-insensitive)
# - lake_col: name of the lake identifier column
exclude_lakes <- function(df, lakes = default_lake_exclusions()$LakeID, lake_col = "LakeID", .warn = TRUE) {
  if (!lake_col %in% names(df)) return(df)

  # Name-based exclusion
  extra_env <- Sys.getenv("EXCLUDE_LAKES_EXTRA", "")
  extra_opt <- getOption("exclude_lakes_extra", NULL)
  extra_vec <- character()
  if (nzchar(extra_env)) extra_vec <- c(extra_vec, strsplit(extra_env, ",")[[1]])
  if (!is.null(extra_opt)) extra_vec <- c(extra_vec, extra_opt)
  extra_vec <- tolower(trimws(extra_vec))
  lakes_all <- unique(c(tolower(lakes), extra_vec))
  lake_vals <- as.character(df[[lake_col]])
  remove_by_name <- tolower(trimws(lake_vals)) %in% lakes_all

  keep <- !remove_by_name
  out <- df[keep, , drop = FALSE]

  if (isTRUE(.warn)) {
    removed <- sum(!keep, na.rm = TRUE)
    removed_names <- unique(tolower(trimws(lake_vals[remove_by_name])))
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
# df <- readRDS("data/main_zooplankton_data.rds")
# df_clean <- exclude_lakes(df)                 # default excludes nimpkish and msk
# ex_tbl <- default_lake_exclusions()           # names + reasons for reporting
# df_clean2 <- exclude_lakes(df, lakes = c("nimpkish","another"))
# # proceed with modeling using df_clean*


