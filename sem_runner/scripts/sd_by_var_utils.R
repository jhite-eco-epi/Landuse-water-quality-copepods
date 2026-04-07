#!/usr/bin/env Rscript
#
# Utilities for computing and using global SDs for standardized (beta*) effects.
#
# Standardization convention used throughout this repo:
#   effect_std = effect_raw * SD(X) / SD(Y)
#
# This module centralizes:
# - computing SDs from a single "canonical" data frame (preferred)
# - a fallback that can compute SDs from multiple model.frames and pick max-n
# - reading/writing sd_by_var.csv in an output directory

sd_by_var_compute <- function(df) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) {
    return(data.frame(var = character(), sd = numeric(), n = integer(), stringsAsFactors = FALSE))
  }
  is_num <- vapply(df, is.numeric, logical(1))
  vars <- names(df)[is_num]
  rows <- vector("list", length(vars))
  k <- 0L
  for (v in vars) {
    x <- suppressWarnings(as.numeric(df[[v]]))
    ok <- is.finite(x)
    n_ok <- sum(ok)
    if (n_ok < 3) next
    s <- stats::sd(x[ok])
    if (!is.finite(s) || s <= 0) next
    k <- k + 1L
    rows[[k]] <- data.frame(var = v, sd = s, n = n_ok, stringsAsFactors = FALSE)
  }
  if (k == 0L) return(data.frame(var = character(), sd = numeric(), n = integer(), stringsAsFactors = FALSE))
  out <- do.call(rbind, rows[seq_len(k)])
  out <- out[order(out$var), , drop = FALSE]
  out
}

sd_by_var_compute_from_frames <- function(frames) {
  # frames: list of data.frames (e.g., model.frame(fit) from multiple submodels)
  if (is.null(frames) || !length(frames)) {
    return(data.frame(var = character(), sd = numeric(), n = integer(), stringsAsFactors = FALSE))
  }
  rows <- list()
  for (f in frames) {
    if (is.null(f) || !is.data.frame(f) || !nrow(f)) next
    rows[[length(rows) + 1L]] <- sd_by_var_compute(f)
  }
  if (!length(rows)) return(data.frame(var = character(), sd = numeric(), n = integer(), stringsAsFactors = FALSE))
  allr <- do.call(rbind, rows)
  if (is.null(allr) || !nrow(allr)) return(data.frame(var = character(), sd = numeric(), n = integer(), stringsAsFactors = FALSE))
  # For each var, pick the SD computed from the largest n (ties: first)
  spl <- split(allr, allr$var)
  picked <- lapply(spl, function(d) d[which.max(d$n), , drop = FALSE])
  out <- do.call(rbind, picked)
  out <- out[order(out$var), , drop = FALSE]
  rownames(out) <- NULL
  out
}

sd_by_var_write <- function(df, out_dir, filename = "sd_by_var.csv") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  tbl <- sd_by_var_compute(df)
  path <- file.path(out_dir, filename)
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_csv(tbl, path)
  } else {
    utils::write.csv(tbl, path, row.names = FALSE)
  }
  invisible(tbl)
}

sd_by_var_read <- function(out_dir, filename = "sd_by_var.csv") {
  path <- file.path(out_dir, filename)
  if (!file.exists(path)) return(NULL)
  if (requireNamespace("readr", quietly = TRUE)) {
    df <- tryCatch(suppressMessages(readr::read_csv(path, show_col_types = FALSE)), error = function(e) NULL)
  } else {
    df <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  }
  if (is.null(df) || !is.data.frame(df) || !nrow(df) || !all(c("var", "sd") %in% names(df))) return(NULL)
  df
}

sd_by_var_lookup_factory <- function(sd_tbl, resp_tag_of = NULL) {
  # Returns a function(var) -> sd or NA_real_
  function(v) {
    if (is.null(sd_tbl) || !is.data.frame(sd_tbl) || !nrow(sd_tbl)) return(NA_real_)
    v <- as.character(v)
    hit <- sd_tbl$sd[match(v, sd_tbl$var)]
    if (is.finite(hit)) return(as.numeric(hit))
    if (!is.null(resp_tag_of) && is.function(resp_tag_of)) {
      vt <- resp_tag_of(v)
      hit <- sd_tbl$sd[match(vt, sd_tbl$var)]
      if (is.finite(hit)) return(as.numeric(hit))
    }
    NA_real_
  }
}

