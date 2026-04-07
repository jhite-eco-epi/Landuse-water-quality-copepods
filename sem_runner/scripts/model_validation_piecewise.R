#!/usr/bin/env Rscript

# Piecewise model validation: per-submodel MCMC diagnostics + PPC and PIT plots.
#
# This script expects an output directory containing `fit_<response>.rds` files.
# It writes:
# - mcmc_diagnostics_summary.csv (one row per submodel)
# - ppc_<response>.png (density overlay PPC)
# - ppc_pit_<response>.csv + ppc_pit_hist_<response>.png + ppc_pit_ecdf_<response>.png
# - ppc_pit_ks_<response>.txt (optional)

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(ggplot2)
  library(posterior)
  library(readr)
})

OUT_DIR <- Sys.getenv("SEM_OUT_DIR", unset = NA_character_)
if (is.na(OUT_DIR) || !nzchar(OUT_DIR) || !dir.exists(OUT_DIR)) {
  stop("SEM_OUT_DIR must point to an existing model output directory. Got: ", OUT_DIR)
}

NDRAWS_PIT <- suppressWarnings(as.integer(Sys.getenv("SEM_VALIDATION_PIT_DRAWS", unset = "200")))
if (!is.finite(NDRAWS_PIT) || NDRAWS_PIT < 50) NDRAWS_PIT <- 200L
NDRAWS_PPC <- suppressWarnings(as.integer(Sys.getenv("SEM_VALIDATION_PPC_DRAWS", unset = "100")))
if (!is.finite(NDRAWS_PPC) || NDRAWS_PPC < 20) NDRAWS_PPC <- 100L

fit_files <- list.files(OUT_DIR, pattern = "^fit_.*\\.rds$", full.names = TRUE)
if (!length(fit_files)) {
  message("No fit files found in ", OUT_DIR, " (expected fit_*.rds). Skipping validation.")
  quit(save = "no", status = 0)
}

safe_int <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return(NA_integer_)
  as.integer(x)
}

safe_num <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return(NA_real_)
  as.numeric(x)
}

diag_rows <- list()

for (fp in sort(fit_files)) {
  fit <- readRDS(fp)
  # For piecewise runs we always save as fit_<response>.rds, so parse response from filename.
  # This avoids brms multiformula edge cases where `formula(fit)` is not a single response.
  resp_var <- sub("^fit_", "", sub("\\.rds$", "", basename(fp)))
  message("Validating: ", resp_var)

  mf <- try(as.data.frame(model.frame(fit)), silent = TRUE)
  n_obs <- if (inherits(mf, "try-error")) NA_integer_ else nrow(mf)

  # --- MCMC diagnostics (R-hat + ESS) ---
  draws <- try(posterior::as_draws_array(fit), silent = TRUE)
  summ <- if (inherits(draws, "try-error")) NULL else {
    posterior::summarise_draws(
      draws,
      rhat = posterior::rhat,
      ess_bulk = posterior::ess_bulk,
      ess_tail = posterior::ess_tail
    )
  }

  max_rhat <- if (is.null(summ) || !("rhat" %in% names(summ))) NA_real_ else suppressWarnings(max(summ$rhat, na.rm = TRUE))
  min_ess_bulk <- if (is.null(summ) || !("ess_bulk" %in% names(summ))) NA_real_ else suppressWarnings(min(summ$ess_bulk, na.rm = TRUE))
  min_ess_tail <- if (is.null(summ) || !("ess_tail" %in% names(summ))) NA_real_ else suppressWarnings(min(summ$ess_tail, na.rm = TRUE))
  n_params <- if (is.null(summ)) NA_integer_ else nrow(summ)

  if (!is.null(summ) && nrow(summ)) {
    out_param <- file.path(OUT_DIR, paste0("mcmc_diagnostics_params_", resp_var, ".csv"))
    readr::write_csv(summ, out_param)
  }

  # --- NUTS diagnostics (divergences / treedepth) ---
  np <- try(brms::nuts_params(fit), silent = TRUE)
  n_divergent <- NA_integer_
  n_treedepth_saturated <- NA_integer_
  if (!inherits(np, "try-error") && is.data.frame(np)) {
    # `nuts_params()` returns long format with columns like: Chain, Iteration, Parameter, Value
    if (all(c("Parameter", "Value") %in% names(np))) {
      n_divergent <- safe_int(sum(np$Parameter == "divergent__" & np$Value == 1, na.rm = TRUE))
      n_treedepth_saturated <- safe_int(sum(np$Parameter == "treedepth__" & np$Value >= 15, na.rm = TRUE))
    }
  }

  # --- Posterior predictive draws (used for both PPC density overlay and PIT) ---
  ppc_path <- file.path(OUT_DIR, paste0("ppc_", resp_var, ".png"))
  pit_csv <- file.path(OUT_DIR, paste0("ppc_pit_", resp_var, ".csv"))
  pit_hist <- file.path(OUT_DIR, paste0("ppc_pit_hist_", resp_var, ".png"))
  pit_ecdf <- file.path(OUT_DIR, paste0("ppc_pit_ecdf_", resp_var, ".png"))
  pit_ks <- file.path(OUT_DIR, paste0("ppc_pit_ks_", resp_var, ".txt"))

  if (!inherits(mf, "try-error") && (resp_var %in% names(mf))) {
    y <- mf[[resp_var]]
    yrep <- try(brms::posterior_predict(fit, ndraws = NDRAWS_PIT), silent = TRUE)
    if (!inherits(yrep, "try-error")) {
      # Density overlay PPC (reuse yrep; no titles)
      p_ppc <- NULL
      if (requireNamespace("bayesplot", quietly = TRUE)) {
        yrep_ppc <- yrep[seq_len(min(NDRAWS_PPC, nrow(yrep))), , drop = FALSE]
        p_ppc <- try(bayesplot::ppc_dens_overlay(y, yrep_ppc), silent = TRUE)
      }
      if (!inherits(p_ppc, "try-error") && !is.null(p_ppc)) {
        p_ppc <- p_ppc + labs(title = NULL, x = resp_var, y = "Density")
        ggsave(ppc_path, p_ppc, width = 7, height = 5, dpi = 150, bg = "white")
      }

      pit <- vapply(seq_len(ncol(yrep)), function(i) mean(yrep[, i] <= y[i]), numeric(1))
      pit_df <- tibble::tibble(idx = seq_along(pit), pit = pit)
      readr::write_csv(pit_df, pit_csv)

      p_hist <- ggplot(pit_df, aes(x = pit)) +
        geom_histogram(bins = 10, fill = "steelblue", color = "white") +
        geom_hline(yintercept = nrow(pit_df) / 10, linetype = "dashed", color = "gray40") +
        coord_cartesian(xlim = c(0, 1)) +
        labs(title = NULL, x = "PIT", y = "Count") +
        theme_minimal()
      ggsave(pit_hist, p_hist, width = 7, height = 5, dpi = 150, bg = "white")

      p_ecdf <- ggplot(pit_df, aes(x = pit)) +
        stat_ecdf(geom = "step", color = "steelblue") +
        stat_function(fun = function(x) x, linetype = "dashed", color = "gray40") +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
        labs(title = NULL, x = "PIT", y = "ECDF") +
        theme_minimal()
      ggsave(pit_ecdf, p_ecdf, width = 7, height = 5, dpi = 150, bg = "white")

      ks <- try(stats::ks.test(pit_df$pit, "punif"), silent = TRUE)
      if (!inherits(ks, "try-error")) {
        capture.output(print(ks), file = pit_ks)
      }
    }
  }

  diag_rows[[length(diag_rows) + 1]] <- tibble::tibble(
    response = resp_var,
    n_obs = safe_int(n_obs),
    n_params = safe_int(n_params),
    max_rhat = safe_num(max_rhat),
    min_ess_bulk = safe_num(min_ess_bulk),
    min_ess_tail = safe_num(min_ess_tail),
    n_divergent = safe_int(n_divergent),
    n_treedepth_saturated = safe_int(n_treedepth_saturated),
    ppc_density = basename(ppc_path),
    pit_csv = basename(pit_csv),
    pit_hist = basename(pit_hist),
    pit_ecdf = basename(pit_ecdf)
  )
}

diag_tbl <- dplyr::bind_rows(diag_rows) %>% dplyr::arrange(.data$response)
readr::write_csv(diag_tbl, file.path(OUT_DIR, "mcmc_diagnostics_summary.csv"))

message("Wrote model validation artifacts to ", OUT_DIR)

