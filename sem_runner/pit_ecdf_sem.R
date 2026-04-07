#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# Paths ------------------------------------------------------------------------
fit_path <- file.path("sem_runner", "results_deep", "sem_gaussian_log1p_rescor.rds")
out_dir  <- file.path("sem_runner", "results_deep")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(fit_path)) stop("Fit not found: ", fit_path)
fit <- readRDS(fit_path)

dat <- fit$data

# Internal brms response tags (underscores removed)
resp_info <- list(
  list(resp = "chllog1p", col = "chl_log1p",  label = "log chlorophyll"),
  list(resp = "cyclog1p", col = "cyc_log1p",  label = "log cyclopoids")
)

# Helper to compute exact LOO-PIT for one response --------------------------------
compute_pit <- function(fit, dat, resp_tag, col_name, ndraws = 1000L) {
  if (!(col_name %in% names(dat))) stop("Column not found in fit$data: ", col_name)
  y <- dat[[col_name]]
  n <- nrow(dat)
  pit <- numeric(n)
  set.seed(123)
  message(sprintf("Computing exact LOO-PIT for %s (n=%d, draws=%d)", resp_tag, n, ndraws))
  for (i in seq_len(n)) {
    train <- dat[-i, , drop = FALSE]
    test  <- dat[i, , drop = FALSE]
    # Refit on training data only; reuse compiled model when possible
    fit_i <- try(
      update(
        fit,
        newdata   = train,
        recompile = TRUE,
        refresh   = 0
      ),
      silent = TRUE
    )
    if (inherits(fit_i, "try-error")) {
      stop("Failed to update/refit SEM model on training data for observation ", i)
    }
    yrep_i <- posterior_predict(fit_i, newdata = test, ndraws = ndraws, resp = resp_tag)
    pit[i] <- mean(as.numeric(yrep_i) <= as.numeric(y[i]))
  }
  tibble(idx = seq_len(n), pit = pit)
}

# Run for both responses --------------------------------------------------------
for (ri in resp_info) {
  resp_tag <- ri$resp
  col_name <- ri$col
  label    <- ri$label
  pit_df <- compute_pit(fit, dat, resp_tag, col_name, ndraws = 1000L)

  # Save CSV
  out_csv <- file.path(out_dir, paste0("loo_pit_exact_", col_name, ".csv"))
  readr::write_csv(pit_df, out_csv)

  # Histogram (no title)
  p_hist <- ggplot(pit_df, aes(x = pit)) +
    geom_histogram(bins = 10, fill = "steelblue", color = "white") +
    geom_hline(yintercept = nrow(pit_df) / 10, linetype = "dashed", color = "gray40") +
    labs(title = NULL, x = "PIT", y = "Count") +
    theme_minimal()
  ggsave(file.path(out_dir, paste0("loo_pit_exact_hist_", col_name, ".png")), p_hist,
         width = 7, height = 5, dpi = 150, bg = "white")

  # ECDF (no title)
  p_ecdf <- ggplot(pit_df, aes(x = pit)) +
    stat_ecdf(geom = "step", color = "steelblue") +
    stat_function(fun = function(x) x, linetype = "dashed", color = "gray40") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = NULL, x = "PIT", y = "ECDF") +
    theme_minimal()
  ggsave(file.path(out_dir, paste0("loo_pit_exact_ecdf_", col_name, ".png")), p_ecdf,
         width = 7, height = 5, dpi = 150, bg = "white")

  # KS test summary
  ks <- try(stats::ks.test(pit_df$pit, "punif"), silent = TRUE)
  if (!inherits(ks, "try-error")) {
    capture.output(print(ks), file = file.path(out_dir, paste0("loo_pit_exact_ks_", col_name, ".txt")))
  }

  message("Saved PIT artifacts for ", label, " to ", out_dir)
}
