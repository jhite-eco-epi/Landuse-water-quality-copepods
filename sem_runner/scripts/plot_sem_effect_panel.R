#!/usr/bin/env Rscript
# Exports a function `make_effect_plot(effect_spec)` which returns a ggplot
# effect plot for a predictor->response pair using the fitted SEM model.
#
# effect_spec format: "predictor->response", e.g. "chl_log1p->cyc_log1p"

suppressPackageStartupMessages({
  library(ggplot2)
  library(brms)
  library(dplyr)
  library(readr)
})

OUT_DIR <- Sys.getenv("SEM_OUT_DIR", unset = "sem_runner/results_deep")
FIT_ENV  <- Sys.getenv("SEM_FIT_RDS", unset = "")
deep_fit_path <- if (nzchar(FIT_ENV)) FIT_ENV else file.path(OUT_DIR, "sem_gaussian_log1p_rescor.rds")
MODEL_NAME <- Sys.getenv("SEM_MODEL_NAME", unset = "gaussian_log1p_rescor")
piecewise_fit_paths <- list.files(OUT_DIR, pattern = "^fit_.*\\.rds$", full.names = TRUE)
use_deep <- file.exists(deep_fit_path)
use_piecewise <- length(piecewise_fit_paths) > 0 && !use_deep
if (!use_deep && !use_piecewise) {
  warning("No SEM fit found in OUT_DIR; effect panels will be placeholders.")
}

friendly <- function(v) {
  switch(v,
    "chl_log1p" = "chlorophyll",
    "cal_log1p" = "calanoids",
    "cyc_log1p" = "cyclopoids",
    "disease"   = "copepod infection burden",
    "pH"        = "pH",
    "oxy"       = "oxygen",
    "temp"      = "temperature",
    "defor"     = "deforestation (% forest loss)",
    "thermo"    = "thermocline depth",
    "MaxDepth"  = "max depth",
    v
  )
}

make_effect_plot <- function(effect_spec) {
  parts <- strsplit(effect_spec, "->", fixed = TRUE)[[1]]
  if (length(parts) != 2) {
    warning("Invalid effect_spec: ", effect_spec)
    return(ggplot() + theme_void())
  }
  predictor <- parts[1]
  response  <- parts[2]
  fit <- NULL
  resp_tag <- switch(response,
    "chl_log1p"  = "chllog1p",
    "cyc_log1p"  = "cyclog1p",
    "cal_log1p"  = "callog1p",
    response
  )
  if (use_deep) {
    fit <- readRDS(deep_fit_path)
  } else if (use_piecewise) {
    fp <- file.path(OUT_DIR, paste0("fit_", response, ".rds"))
    if (file.exists(fp)) fit <- readRDS(fp)
  }
  if (is.null(fit)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "SEM fit not found", size = 6) + theme_void())
  }

  # Prefer conditional_effects because it works on the response scale for non-Gaussian families
  ce <- try({
    if (use_deep) conditional_effects(fit, resp = resp_tag, effects = predictor) else conditional_effects(fit, effects = predictor)
  }, silent = TRUE)
  if (inherits(ce, "try-error")) {
    # Fallback: use coefficients CSV (linear/identity approximation)
    coef_path <- file.path(OUT_DIR, sprintf("coefficients_%s_%s.csv", resp_tag, MODEL_NAME))
    if (file.exists(coef_path)) {
      tbl <- tryCatch(readr::read_csv(coef_path, show_col_types = FALSE), error = function(e) NULL)
      if (!is.null(tbl)) {
        intercept_row <- tbl[grepl("_Intercept$", tbl$Parameter), , drop = FALSE]
        slope_row     <- tbl[grepl(paste0("_", predictor, "$"), tbl$Parameter), , drop = FALSE]
        if (nrow(slope_row) == 1 && nrow(intercept_row) == 1) {
          x_obs <- suppressWarnings(as.numeric(fit$data[[predictor]]))
          x_rng <- if (all(is.finite(x_obs))) stats::quantile(x_obs, probs = c(0.02, 0.98), na.rm = TRUE) else c(-2, 2)
          x_seq <- seq(x_rng[1], x_rng[2], length.out = 200)
          y_mean  <- intercept_row$Mean + slope_row$Mean * x_seq
          y_lower <- intercept_row$Q2.5 + slope_row$Q2.5 * x_seq
          y_upper <- intercept_row$Q97.5 + slope_row$Q97.5 * x_seq
          df <- tibble::tibble(x = x_seq, y = y_mean, ylo = y_lower, yhi = y_upper)
          gp <- ggplot(df, aes(x = x, y = y)) +
            geom_ribbon(aes(ymin = ylo, ymax = yhi), fill = "grey70", alpha = 0.3, color = NA) +
            geom_line(color = "#2C7FB8", linewidth = 1) +
            labs(title = NULL, x = friendly(predictor), y = friendly(response)) +
            theme_classic() +
            theme(
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
              panel.background = element_blank(),
              plot.background = element_blank(),
              panel.grid = element_blank(),
              axis.line = element_blank()
            )
          return(gp)
        }
      }
    }
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste0("No effect: ", effect_spec), size = 5) + theme_void())
  }
  pls <- try(plot(ce, plot = FALSE), silent = TRUE)
  if (inherits(pls, "try-error") || length(pls) == 0) {
    return(ggplot() + theme_void())
  }
  gp <- pls[[1]] +
    # Add observed points (lake-level medians used in the SEM fit) for context.
    # This avoids relying on brms' plot defaults and keeps the effect panels consistent.
    {
      if (all(c(predictor, response) %in% names(fit$data))) {
        obs <- fit$data %>%
          dplyr::select(dplyr::all_of(c(predictor, response))) %>%
          dplyr::rename(x = dplyr::all_of(predictor), y = dplyr::all_of(response)) %>%
          dplyr::filter(is.finite(.data$x), is.finite(.data$y))
        ggplot2::geom_point(
          data = obs,
          mapping = ggplot2::aes(x = .data$x, y = .data$y),
          inherit.aes = FALSE,
          alpha = 0.35,
          size = 1.6,
          color = "black"
        )
      } else {
        NULL
      }
    } +
    labs(title = NULL, x = friendly(predictor), y = friendly(response)) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_blank()
    )
  return(gp)
}


