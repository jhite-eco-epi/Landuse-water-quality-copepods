#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

OUT_DIR  <- Sys.getenv("SEM_OUT_DIR", unset = "sem_runner/results_deep")
MODEL    <- Sys.getenv("SEM_MODEL_NAME", unset = "model_run")
FIT_PATH <- Sys.getenv("SEM_FIT_RDS", unset = "")
SPEC_ENV <- Sys.getenv("SEM_SPEC_PATH", unset = "")

if (!dir.exists(OUT_DIR)) stop("Output directory not found: ", OUT_DIR)

md <- c()
emit <- function(...) { md <<- c(md, paste0(...)) }
emitln <- function(...) { emit(paste0(..., "\n")) }

emitln("# SEM model: ", MODEL)
emitln()
emitln("Output directory: `", OUT_DIR, "`")
if (nzchar(FIT_PATH)) emitln("Fit file: `", basename(FIT_PATH), "`")
emitln()

fmt_int <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return("?")
  format(round(as.numeric(x)), big.mark = ",", trim = TRUE)
}

# Helper: markdown table from data frame (simple)
md_table <- function(df, max_rows = 10) {
  df2 <- head(df, max_rows)
  cols <- names(df2)
  header <- paste0("| ", paste(cols, collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", length(cols)), collapse = "|"), "|")
  rows <- apply(df2, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  c(header, sep, rows)
}

ctx <- NULL
ctx_path <- file.path(OUT_DIR, "readme_context.rds")
if (file.exists(ctx_path)) {
  ctx <- try(readRDS(ctx_path), silent = TRUE)
  if (inherits(ctx, "try-error")) ctx <- NULL
}

if (!is.null(ctx)) {
  emitln("## Lake coverage")
  emitln()
  if (!is.null(ctx$final_lakes) && !is.na(ctx$final_lakes)) {
    emitln("Model uses ", fmt_int(ctx$final_lakes), " lake-level medians (deep/single habitat).")
  } else if (!is.null(ctx$piecewise_n_by_equation) && is.data.frame(ctx$piecewise_n_by_equation) && nrow(ctx$piecewise_n_by_equation)) {
    pwt0 <- ctx$piecewise_n_by_equation
    rng <- range(pwt0$n_lakes, na.rm = TRUE)
    if (all(is.finite(rng))) {
      emitln("Piecewise fits use varying sample sizes across equations: n = ", fmt_int(rng[1]), "–", fmt_int(rng[2]), " lakes.")
    }
  }
  if (!is.null(ctx$deep_single_candidates)) {
    emitln("Deep/single candidates prior to dropping missing covariates: ", fmt_int(ctx$deep_single_candidates), ".")
  }
  if (!is.null(ctx$total_lakes_raw)) {
    emitln("Unique lakes in the main dataset before exclusions: ", fmt_int(ctx$total_lakes_raw), ".")
  }
  if (!is.null(ctx$max_depth_threshold)) {
    emitln("Rule: lakes with MaxDepth <", ctx$max_depth_threshold, " m are excluded.")
  }
  emitln()

  # Optional: piecewise SEM per-equation sample sizes (when fits use different n)
  if (!is.null(ctx$piecewise_n_by_equation) && is.data.frame(ctx$piecewise_n_by_equation) && nrow(ctx$piecewise_n_by_equation)) {
    emitln("### Piecewise per-equation sample sizes")
    emitln()
    pwt <- ctx$piecewise_n_by_equation
    # Keep a compact view for README; full table is written as CSV in the output directory.
    pwt <- pwt[order(-pwt$n_lakes), , drop = FALSE]
    if (nrow(pwt) > 15) pwt <- pwt[1:15, , drop = FALSE]
    emitln(paste(md_table(pwt, max_rows = 15), collapse = "\n"))
    emitln()
  }

  # Optional: piecewise missing-data analysis per submodel
  if (!is.null(ctx$piecewise_missingness_summary) && is.data.frame(ctx$piecewise_missingness_summary) && nrow(ctx$piecewise_missingness_summary)) {
    emitln("### Missing data by submodel (piecewise)")
    emitln()
    emitln("Each submodel is fit on the subset of lakes with complete data for the response and its predictors. The table below summarizes sample sizes and the most common missing variables; per-submodel dropped-lake lists are written as `missingness_<response>.csv` in the output directory.")
    emitln()
    ms <- ctx$piecewise_missingness_summary
    ms <- ms %>%
      dplyr::select(equation, response, n_total, n_used, n_dropped, top_missing) %>%
      dplyr::arrange(desc(.data$n_dropped), .data$equation)
    emitln(paste(md_table(ms, max_rows = 50), collapse = "\n"))
    emitln()
  }
}

# Model validation (piecewise)
diag_path <- file.path(OUT_DIR, "mcmc_diagnostics_summary.csv")
if (file.exists(diag_path)) {
  emitln("## Model validation (MCMC + posterior predictive checks)")
  emitln()
  emitln("For each piecewise submodel we record basic MCMC convergence diagnostics (R-hat and effective sample sizes) and run posterior predictive checks (density overlays) plus posterior predictive PIT checks (histogram + ECDF vs Uniform(0,1)). PIT here is computed from posterior predictive draws on the fitted data (not leave-one-out).")
  emitln()
  diag <- tryCatch(suppressMessages(readr::read_csv(diag_path, show_col_types = FALSE)), error = function(e) NULL)
  if (!is.null(diag) && nrow(diag)) {
    keep <- intersect(names(diag), c("response","n_obs","max_rhat","min_ess_bulk","min_ess_tail","n_divergent","n_treedepth_saturated"))
    if (length(keep)) diag <- diag[, keep, drop = FALSE]
    emitln(paste(md_table(diag, max_rows = 100), collapse = "\n"))
    emitln()
    emitln("Artifacts (per response) are written to the output directory with prefixes `ppc_` (density PPC) and `ppc_pit_` (PIT histogram/ECDF), and a full parameter-level diagnostics table `mcmc_diagnostics_params_<response>.csv`.")
    emitln()
  }
}

# Resolve composite effects list (for legend)
EFFECT_LIST <- character()
if (nzchar(SPEC_ENV) && file.exists(SPEC_ENV)) {
  try(sys.source(SPEC_ENV, envir = environment()), silent = TRUE)
}
if (exists("EFFECT_PANELS", inherits = FALSE) && is.character(EFFECT_PANELS)) {
  EFFECT_LIST <- EFFECT_PANELS
} else if (exists("PATHS", inherits = FALSE) && length(PATHS) > 0) {
  term_of <- function(node) {
    switch(node,
      "chllog1p" = "chl_log1p",
      "callog1p" = "cal_log1p",
      "cyclog1p" = "cyc_log1p",
      node
    )
  }
  dl <- character()
  for (pth in PATHS) {
    if (length(pth) < 2) next
    for (k in seq_len(length(pth) - 1)) {
      dl <- c(dl, paste0(term_of(pth[k]), "->", term_of(pth[k + 1])))
    }
  }
  EFFECT_LIST <- unique(dl)
}
friendly <- function(v) {
  switch(v,
    "chl_log1p" = "chlorophyll",
    "cal_log1p" = "calanoids",
    "cyc_log1p" = "cyclopoids",
    "eutroph"   = "eutrophication composite",
    "thermo"    = "thermocline depth",
    "MaxDepth"  = "max depth",
    "pH"        = "pH",
    "oxy"       = "oxygen",
    "temp"      = "temperature",
    "defor"     = "deforestation (% forest loss)",
    v
  )
}

preprocessing_exclusion_text <- function(ctx) {
  depth_txt <- ""
  if (!is.null(ctx$max_depth_threshold)) {
    depth_txt <- paste0(" Lakes with MaxDepth < ", ctx$max_depth_threshold, " m were excluded.")
  }
  depth_txt
}

# Composite figure and legend
comp_png <- file.path(OUT_DIR, "composite_sem_AA_BC.png")
if (file.exists(comp_png)) {
  path_nodes <- character()
  if (exists("PATHS", inherits = FALSE) && length(PATHS) > 0) {
    path_nodes <- unique(unlist(PATHS))
  }
  node_lines <- c(
    "thermocline depth = “thermo”",
    "chlorophyll = “chl_log1p”",
    "calanoids = “cal_log1p”",
    "cyclopoids = “cyc_log1p”",
    "pH = “pH”",
    "oxygen = “oxy”"
  )
  if ("temp" %in% path_nodes) node_lines <- c(node_lines, "temperature = “temp”")
  if ("defor" %in% path_nodes) node_lines <- c(node_lines, "deforestation = “defor”")
  if ("eutroph" %in% path_nodes) node_lines <- c(node_lines, "eutrophication composite = “eutroph”")
  emitln("## Composite figure")
  emitln()
  emitln("![](", basename(comp_png), ")")
  emitln()
  # Paragraph‑style, journal‑ready legend
  rescor_txt <- if (exists("rescor_flag", inherits = FALSE)) {
    if (isTRUE(rescor_flag)) "Residual correlations were estimated where available." else "Residual correlations were set to zero (rescor = FALSE)."
  } else {
    "Residual correlations were set to zero (rescor = FALSE)."
  }
  emitln("Figure. Piecewise structural equation model (SEM) of deep/single lakes. The path diagram in panel A shows directed links among responses; edge labels are posterior median regression coefficients from BRMS Gaussian components (density variables on the natural-log scale, ln(1 + x)). ", rescor_txt)
  emitln()
  emitln("Panel A (Path diagram). Directed graph with edges labeled by posterior median standardized coefficients (global SD scaling). Raw and standardized edge summaries are saved in `path_diagram_edges.csv` (columns `*_raw` and `*_std`). Node names: ", paste(node_lines, collapse = ", "), ".")
  if (length(EFFECT_LIST) > 0) {
    last_letter <- LETTERS[1 + length(EFFECT_LIST)]
    # Build a compact enumeration “(B) x→y; (C) …”
    parts_txt <- c()
    for (i in seq_along(EFFECT_LIST)) {
      pr <- strsplit(EFFECT_LIST[[i]], "->", fixed = TRUE)[[1]]
      if (length(pr) == 2) {
        parts_txt <- c(parts_txt, paste0("(", LETTERS[i + 1], ") ", friendly(pr[1]), " → ", friendly(pr[2])) )
      }
    }
    if (length(parts_txt)) {
      emitln()
      emitln("Panels B–", last_letter, " (conditional effects). Each panel displays the fitted mean response (solid line) as a function of the listed predictor, with a 95% credible band (2.5%–97.5% posterior quantiles of the intercept and slope) shown in gray; x‑axes span the 2nd–98th percentiles of the observed predictor. Paths shown: ", paste(parts_txt, collapse = "; "), ".")
    }
  } else {
    emitln()
    emitln("Panels B–C (conditional effects). Fitted mean responses with 95% credible bands for exemplar paths; x‑axes cover the central 96% of observed predictor values (2nd–98th percentiles).")
  }
  emitln()
  emitln("Data and preprocessing. All abiotic predictors are water‑column integrated summaries (integrated_chl, integrated_fDOM, integrated_pH, integrated_temp, integrated_DO_percent; thermocline_depth_m from profiles). Density responses are modeled on the natural-log scale using ln(1 + x) transforms. ", preprocessing_exclusion_text(ctx))
  emitln()
}

# Coefficient tables
coef_files <- list.files(OUT_DIR, pattern = "^coefficients_.*\\.csv$", full.names = TRUE)
if (length(coef_files)) {
  emitln("## Coefficients")
  for (cf in sort(coef_files)) {
    emitln()
    emitln("### ", basename(cf))
    df <- tryCatch(suppressMessages(readr::read_csv(cf, show_col_types = FALSE)), error = function(e) NULL)
    if (!is.null(df) && nrow(df)) {
      emitln(paste(md_table(df), collapse = "\n"))
      emitln()
    } else {
      emitln("_Empty or unreadable_: `", basename(cf), "`")
    }
  }
  emitln()
}

# Path effects (direct / indirect / total)
pf_direct <- file.path(OUT_DIR, "path_effects.csv")
pf_ind    <- file.path(OUT_DIR, "path_effects_indirect.csv")
pf_total  <- file.path(OUT_DIR, "path_effects_total.csv")
if (file.exists(pf_direct) || file.exists(pf_total) || file.exists(pf_ind)) {
  emitln("## Effects (direct, indirect, total)")
  emitln()
  emitln("Effects are computed from posterior draws. Direct effects are regression slopes for each arrow. Indirect effects are computed draw-by-draw as products of coefficients along each directed path; total effects are direct + summed indirect (when applicable).")
  emitln()

  if (file.exists(pf_total)) {
    emitln("### Total effects (`path_effects_total.csv`)")
    df <- tryCatch(suppressMessages(readr::read_csv(pf_total, show_col_types = FALSE)), error = function(e) NULL)
    if (!is.null(df) && nrow(df)) {
      # Keep the README compact
      keep <- intersect(names(df), c("from_label","to_label","direct_included","n_indirect_paths","nonzero_95","median","q2.5","q97.5","median_std","q2.5_std","q97.5_std","nonzero_95_std","n_draws"))
      if (length(keep)) df <- df[, keep, drop = FALSE]
      emitln(paste(md_table(df, max_rows = 60), collapse = "\n"))
      emitln()
    }
  }

  if (file.exists(pf_direct)) {
    emitln("### Direct effects (`path_effects.csv`)")
    df <- tryCatch(suppressMessages(readr::read_csv(pf_direct, show_col_types = FALSE)), error = function(e) NULL)
    if (!is.null(df) && nrow(df)) {
      # keep compact columns if standardized outputs exist
      keep <- intersect(names(df), c("path","median","q2.5","q97.5","nonzero_95","median_std","q2.5_std","q97.5_std","nonzero_95_std","n_draws"))
      if (length(keep)) df <- df[, keep, drop = FALSE]
      emitln(paste(md_table(df, max_rows = 60), collapse = "\n"))
      emitln()
    }
  }

  if (file.exists(pf_ind)) {
    emitln("### Indirect effects (`path_effects_indirect.csv`)")
    emitln()
    emitln("Indirect effects are reported per directed path (with a `via` and explicit `path`). See `", basename(pf_ind), "` for the full table.")
    emitln()
  }

  emitln()
}

# Figures
fig_groups <- list(
  "Posterior betas" = "^posteriors_betas_.*\\.png$",
  "Posterior sigmas" = "^posteriors_sigma\\.png$",
  "Posterior residual correlation" = "^posteriors_rescor\\.png$",
  "Posterior predictive checks" = "^ppc_.*\\.png$",
  "Path diagram" = "^path_diagram\\.png$"
)
emitln("## Figures")
for (ttl in names(fig_groups)) {
  patt <- fig_groups[[ttl]]
  figs <- list.files(OUT_DIR, pattern = patt, full.names = TRUE)
  if (!length(figs)) next
  emitln()
  emitln("### ", ttl)
  for (fp in sort(figs)) {
    emitln("![](", basename(fp), ")")
  }
}
emitln()

# Write README.md
out_md <- file.path(OUT_DIR, "README.md")
writeLines(md, con = out_md)
message("Wrote ", out_md)


