#!/usr/bin/env Rscript

# Piecewise SEM runner: fit each equation as a separate brms model
# so each "piece" can use all available data (different n per equation).
#

suppressPackageStartupMessages({
  library(brms)
  library(tidyverse)
  library(posterior)
  library(ggplot2)
  library(patchwork)
  library(readr)
})

options(mc.cores = parallel::detectCores())
try({ rstan::rstan_options(auto_write = TRUE, open_progress = TRUE) }, silent = TRUE)

if (identical(Sys.getenv("SEM_DEBUG_TRACEBACK", unset = "0"), "1")) {
  options(error = function() {
    traceback(30)
    quit(save = "no", status = 1)
  })
}

message("=== SEM RUNNER (PIECEWISE BRMS; PER-EQUATION DATA): START ===")

# -------------------- CONFIG --------------------
DATA_RDS <- "data/main_zooplankton_data.rds"
OUT_DIR_BASE <- "sem_runner/results_piecewise"
SEED <- 42
# ------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key) {
  hit <- args[startsWith(args, key)]
  if (length(hit) == 0) return(NA_character_) else sub(paste0("^", key), "", hit[[1]])
}
get_flag <- function(key, default = FALSE) {
  val <- get_arg(key)
  if (is.na(val) || !nzchar(val)) return(default)
  tolower(val) %in% c("1", "true", "t", "yes", "y")
}
spec_path <- get_arg("--spec=")
if (is.na(spec_path) || spec_path == "") {
  spec_env <- Sys.getenv("SEM_SPEC_PATH", unset = NA_character_)
  spec_path <- if (!is.na(spec_env) && nzchar(spec_env)) spec_env else NA_character_
}
keep_all_lakes <- get_flag("--keep-all-lakes=", default = FALSE)
if (is.na(spec_path) || !file.exists(spec_path)) {
  stop("A model spec is required. Pass with --spec=sem_runner/specs/sem_<name>.R or set SEM_SPEC_PATH.")
}

spec_txt <- tryCatch(paste(readLines(spec_path, warn = FALSE), collapse = "\n"), error = function(e) "")
# Determine model_name from spec text without evaluating it
model_name <- NA_character_
m <- regexec("MODEL_NAME\\s*<-\\s*['\\\"]([^'\\\"]+)['\\\"]", spec_txt)
mm <- regmatches(spec_txt, m)[[1]]
if (length(mm) >= 2 && nzchar(mm[2])) model_name <- mm[2]
if (is.na(model_name) || !nzchar(model_name)) {
  bn <- basename(spec_path)
  model_name <- sub("^sem_(.+)\\.R$", "\\1", bn)
}
if (!nzchar(model_name)) stop("Could not determine model name. Define MODEL_NAME in the spec or name file as sem_<model_name>.R.")

OUT_DIR <- file.path(OUT_DIR_BASE, model_name)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(DATA_RDS)) stop("Main dataset not found: ", DATA_RDS)
POSTPROCESSING_SCRIPT <- "data_pipeline/postprocessing.R"
if (!file.exists(POSTPROCESSING_SCRIPT)) stop("Postprocessing helpers not found: ", POSTPROCESSING_SCRIPT)
source(POSTPROCESSING_SCRIPT)
df_raw <- readRDS(DATA_RDS)

safe_num <- function(x) suppressWarnings(as.numeric(x))
has_col <- function(nm) nm %in% names(df_raw)
first_present <- function(cands) {
  nm <- cands[cands %in% names(df_raw)][1]
  if (length(nm) == 0 || is.na(nm)) return(NA_character_) else nm
}
require_col <- function(nm, label) {
  if (!has_col(nm)) stop("Missing required integrated column for ", label, ": ", nm)
  nm
}

resp_tag <- function(v) {
  switch(v,
    "chl_log1p" = "chllog1p",
    "cal_log1p" = "callog1p",
    "cyc_log1p" = "cyclog1p",
    v
  )
}

friendly <- function(v) {
  switch(v,
    "chl_log1p" = "chlorophyll",
    "chllog1p"  = "chlorophyll",
    "cyc_log1p" = "cyclopoids",
    "cyclog1p"  = "cyclopoids",
    "cal_log1p" = "calanoids",
    "callog1p"  = "calanoids",
    "eutroph"   = "eutrophication composite",
    "disease"   = "copepod infection burden",
    "thermo"    = getOption("sem_label_thermo", "thermocline depth"),
    "MaxDepth"  = "max depth",
    "fdom"      = "fDOM",
    "pH"        = "pH",
    "oxy"       = "oxygen",
    "temp"      = "temperature",
    "defor"     = "deforestation (% forest loss)",
    v
  )
}

# -------------------- Preprocess to lake-level medians --------------------
COL_LAKE  <- first_present(c("LakeID", "lake_id", "lakeID"))
COL_HAB   <- first_present(c("zoop_sample_type", "habitat"))
COL_CHL   <- require_col("integrated_chl", "chlorophyll")
COL_PH    <- require_col("integrated_pH", "pH")
COL_CYC   <- first_present(c("Cyclopoids_density_per_L", "cyclopoids", "Cyclopoids"))
COL_THERMO<- first_present(c("thermocline_depth_m", "thermocline_depth", "thermocline_m"))
COL_FDOM  <- require_col("integrated_fDOM", "fDOM")
COL_OXY   <- require_col("integrated_DO_percent", "oxygen (%)")
COL_TEMP  <- first_present(c("integrated_temp", "final_temp", "surface_temp", "bottom_temp", "fallback_temp"))
COL_CAL   <- first_present(c(
  "Calanoids_density_per_L","Calanoid_density_per_L","calanoid_density_per_L",
  "calanoids","Calanoids","Calanoid","calanoid","Calanoida","cal"
))
COL_DEFOR <- first_present(c(
  "mean_crown_cover","mean_crown_cover_pct","crown_cover_mean","crowncover_mean",
  "forest_cover_mean","forest_cover","forest_percent","forest_pct",
  "mean_tree_cover","tree_cover_mean","percent_tree_cover","pct_tree_cover",
  "percent_crown_cover","pct_crown_cover","crown_cover_percent","crown_cover_pct",
  "tree_canopy_percent","tree_canopy_pct","nlcd_tree_canopy","nlcd_forest",
  "deforestation","deforest","crown_cover","landuse_mean_crown_closure"
))

need <- function(nm, label) if (is.na(nm)) stop("Missing required column for ", label)
need(COL_LAKE, "LakeID")
need(COL_HAB, "habitat")
need(COL_CYC, "cyclopoids")
need(COL_THERMO, "thermocline depth")

# exclusions aligned with shared postprocessing defaults
manual_exclusions <- default_lake_exclusions()

df <- df_raw
total_lakes_raw <- length(unique(df[[COL_LAKE]]))
manual_hits <- manual_exclusions[0, , drop = FALSE]
if (!keep_all_lakes && nrow(manual_exclusions)) {
  lake_ids_norm <- tolower(trimws(as.character(df[[COL_LAKE]])))
  manual_hits <- manual_exclusions %>%
    dplyr::mutate(.lake_norm = tolower(trimws(as.character(LakeID)))) %>%
    dplyr::filter(.lake_norm %in% unique(lake_ids_norm)) %>%
    dplyr::select(-.lake_norm)
  df <- exclude_lakes(df, lakes = manual_exclusions$LakeID, lake_col = COL_LAKE, .warn = FALSE)
}

df <- df %>% dplyr::filter(.data[[COL_HAB]] %in% c("deep", "single"))
deep_single_candidates <- length(unique(df[[COL_LAKE]]))

df_num <- df %>%
  transmute(
    LakeID = .data[[COL_LAKE]],
    chl      = safe_num(.data[[COL_CHL]]),
    pH       = safe_num(.data[[COL_PH]]),
    cyc      = safe_num(.data[[COL_CYC]]),
    thermo   = safe_num(.data[[COL_THERMO]]),
    fdom     = safe_num(.data[[COL_FDOM]]),
    oxy      = safe_num(.data[[COL_OXY]]),
    temp     = if (!is.na(COL_TEMP)) safe_num(.data[[COL_TEMP]]) else NA_real_,
    cal      = if (!is.na(COL_CAL)) safe_num(.data[[COL_CAL]]) else NA_real_,
    crown_closure = if (!is.na(COL_DEFOR)) safe_num(.data[[COL_DEFOR]]) else NA_real_
  )

agg <- df_num %>%
  group_by(LakeID) %>%
  summarise(
    chl   = median(chl,   na.rm = TRUE),
    pH    = median(pH,    na.rm = TRUE),
    cyc   = median(cyc,   na.rm = TRUE),
    thermo= median(thermo,na.rm = TRUE),
    fdom  = median(fdom,  na.rm = TRUE),
    oxy   = median(oxy,   na.rm = TRUE),
    temp  = median(temp,  na.rm = TRUE),
    cal   = median(cal,   na.rm = TRUE),
    crown_closure = median(crown_closure, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    defor = {
      fl <- ifelse(is.finite(crown_closure), 100 - crown_closure, NA_real_)
      pmin(100, pmax(0, fl))
    },
    chl_log1p = log1p(pmax(0, chl)),
    cyc_log1p = log1p(pmax(0, cyc)),
    cal_log1p = log1p(pmax(0, cal))
  )

# Allow specs to derive lake-level variables in the active runner environment.
# We alias `agg` as `df_m` so existing deep-runner spec transforms can be reused.
df_m <- agg
sys.source(spec_path, envir = environment())
if (exists("df_m", inherits = FALSE) && is.data.frame(df_m)) {
  agg <- df_m
}

# Global SD table for standardization (computed on deep/single lake medians)
source("sem_runner/scripts/sd_by_var_utils.R", local = TRUE)
sd_by_var_write(agg, OUT_DIR)

core_vars <- c("chl_log1p", "thermo", "fdom", "pH", "oxy", "cyc_log1p")
agg_core <- agg %>% dplyr::filter(if_all(all_of(core_vars), ~ is.finite(.)))
if (nrow(agg_core) < 8) stop("Too few complete deep lakes for core SEM variables.")

# Grab configuration from spec if provided
rescor_flag <- if (exists("rescor_flag", inherits = FALSE)) get("rescor_flag", inherits = FALSE) else FALSE
PATHS <- if (exists("PATHS", inherits = FALSE)) get("PATHS", inherits = FALSE) else list()
EFFECT_PANELS <- if (exists("EFFECT_PANELS", inherits = FALSE)) get("EFFECT_PANELS", inherits = FALSE) else character()
pri_all <- if (exists("pri", inherits = FALSE)) get("pri", inherits = FALSE) else NULL

# Collect piecewise equations: all plain brmsformula objects starting with "bf_"
bf_names <- ls(envir = environment(), pattern = "^bf_", all.names = TRUE)
bf_objs <- lapply(bf_names, function(nm) get(nm, envir = environment()))
names(bf_objs) <- bf_names

is_piece <- function(x) inherits(x, "brmsformula") && !inherits(x, "mvbrmsformula")
bf_piece <- bf_objs[vapply(bf_objs, is_piece, logical(1))]

# Drop combined multivariate holders if present
bf_piece <- bf_piece[setdiff(names(bf_piece), c("bf_chl"))]

if (length(bf_piece) == 0) {
  stop("Spec did not define any piecewise equations (expected brmsformula objects like bf_temp, bf_oxy, ...).")
}

get_response_var <- function(bf) {
  f <- stats::formula(bf)
  as.character(f)[2]
}

ctrl <- list(adapt_delta = 0.99, max_treedepth = 15)

fit_piece <- list()
coverage <- list()
missingness_files <- list()
missingness_summary <- list()

needed_vars_for_bf <- function(bf) {
  f <- stats::formula(bf)
  vars <- unique(all.vars(f))
  # Drop spurious "1" and keep only those present in agg
  vars <- setdiff(vars, "1")
  vars[vars %in% names(agg)]
}

missingness_for_vars <- function(vars) {
  if (!length(vars)) return(data.frame(LakeID = agg$LakeID, missing_vars = "", stringsAsFactors = FALSE))
  miss <- agg %>%
    dplyr::select(LakeID, dplyr::all_of(vars))
  miss_mat <- as.matrix(dplyr::select(miss, dplyr::all_of(vars)))
  miss_labs <- apply(miss_mat, 1, function(row) {
    m <- !is.finite(as.numeric(row))
    if (!any(m)) return("")
    paste(vars[m], collapse = ",")
  })
  data.frame(
    LakeID = miss$LakeID,
    missing_vars = miss_labs,
    stringsAsFactors = FALSE
  )
}

message("Fitting ", length(bf_piece), " piecewise brms models...")
for (nm in names(bf_piece)) {
  message("Preparing equation: ", nm)
  bf <- bf_piece[[nm]]
  resp_var <- get_response_var(bf)
  message("  response var: ", resp_var)
  vars_needed <- needed_vars_for_bf(bf)
  miss_tbl0 <- missingness_for_vars(vars_needed)
  used_ids <- miss_tbl0$LakeID[miss_tbl0$missing_vars == ""]
  dropped <- miss_tbl0 %>%
    dplyr::filter(.data$missing_vars != "") %>%
    dplyr::mutate(
      missing_friendly = vapply(strsplit(.data$missing_vars, ",", fixed = TRUE), function(vs) {
        vs <- trimws(vs)
        paste(vapply(vs, friendly, character(1)), collapse = ", ")
      }, character(1))
    ) %>%
    dplyr::select(LakeID, missing_friendly)

  miss_path <- file.path(OUT_DIR, paste0("missingness_", resp_var, ".csv"))
  readr::write_csv(dropped, miss_path)
  missingness_files[[resp_var]] <- basename(miss_path)

  # Summarize missingness reasons by variable for this equation
  miss_counts <- tibble::tibble(var = vars_needed) %>%
    dplyr::mutate(
      missing_n = vapply(vars_needed, function(v) sum(!is.finite(agg[[v]])), numeric(1)),
      label = vapply(vars_needed, friendly, character(1))
    ) %>%
    dplyr::arrange(desc(.data$missing_n))

  missingness_summary[[length(missingness_summary) + 1]] <- tibble::tibble(
    equation = nm,
    response = resp_var,
    n_total = nrow(agg),
    n_used = length(used_ids),
    n_dropped = nrow(agg) - length(used_ids),
    vars_required = paste(vars_needed, collapse = ", "),
    top_missing = paste0(head(miss_counts$label, 4), " (", head(miss_counts$missing_n, 4), ")", collapse = "; ")
  )
  fit_path <- file.path(OUT_DIR, paste0("fit_", resp_var, ".rds"))
  if (file.exists(fit_path)) {
    message("  cached fit found; loading: ", basename(fit_path))
    fit0 <- readRDS(fit_path)
    fit_piece[[resp_var]] <- fit0
    n0 <- tryCatch(nrow(as.data.frame(model.frame(fit0))), error = function(e) NA_integer_)
    coverage[[length(coverage) + 1]] <- data.frame(
      equation = nm, response = resp_var, n_lakes = n0, stringsAsFactors = FALSE
    )
    next
  }
  # Build minimal model.frame to get complete cases for this equation
  mf <- model.frame(stats::formula(bf), data = agg, na.action = na.omit)
  mf <- as.data.frame(mf)
  message("  model.frame rows: ", nrow(mf))
  n_lakes <- nrow(mf)
  coverage[[length(coverage) + 1]] <- data.frame(
    equation = nm,
    response = resp_var,
    n_lakes = n_lakes,
    stringsAsFactors = FALSE
  )
  message("  coverage row recorded")
  if (n_lakes < 8) {
    warning("Skipping ", nm, " (response ", resp_var, ") due to too few complete lakes: ", n_lakes)
    next
  }

  # Use brms defaults for piecewise fits (safe + avoids resp-scoped priors).
  pri_use <- NULL

  message("- ", nm, " | response=", resp_var, " | n=", n_lakes)

  fit <- brm(
    bf,
    data = mf,
    prior = pri_use,
    chains = 4, iter = 4000, warmup = 1000,
    seed = SEED,
    control = ctrl,
    backend = "rstan",
    refresh = 0,
    save_pars = save_pars(all = TRUE)
  )

  fit_piece[[resp_var]] <- fit
  saveRDS(fit, file.path(OUT_DIR, paste0("fit_", resp_var, ".rds")))
}

coverage_tbl <- dplyr::bind_rows(coverage) %>% arrange(desc(n_lakes))
readr::write_csv(coverage_tbl, file.path(OUT_DIR, "piecewise_n_by_equation.csv"))

missingness_summary_tbl <- dplyr::bind_rows(missingness_summary)
readr::write_csv(missingness_summary_tbl, file.path(OUT_DIR, "piecewise_missingness_summary.csv"))

if (length(fit_piece) == 0) stop("No piecewise fits completed.")

# -------------------- Model validation (PPC + PIT + MCMC diagnostics) --------------------
try({
  val_script <- "sem_runner/scripts/model_validation_piecewise.R"
  if (file.exists(val_script)) {
    Sys.setenv(SEM_OUT_DIR = OUT_DIR)
    message("Sourcing model validation script: ", val_script)
    sys.source(val_script, envir = environment())
  }
}, silent = FALSE)

# -------------------- Coefficient tables (per response) --------------------
coef_summary <- function(x) {
  x <- as.numeric(x)
  tibble(
    Mean = mean(x, na.rm = TRUE),
    SD = stats::sd(x, na.rm = TRUE),
    Q2.5 = as.numeric(stats::quantile(x, 0.025, na.rm = TRUE)),
    Q97.5 = as.numeric(stats::quantile(x, 0.975, na.rm = TRUE))
  )
}

for (resp_var in names(fit_piece)) {
  fit <- fit_piece[[resp_var]]
  dra <- posterior::as_draws_df(fit)
  cols <- grep("^b_", names(dra), value = TRUE)
  if (!length(cols)) next
  rows <- lapply(cols, function(nm) {
    term <- sub("^b_", "", nm)
    out <- coef_summary(dra[[nm]])
    cbind(
      data.frame(Parameter = nm, Term = if (term == "Intercept") "Intercept" else term, stringsAsFactors = FALSE),
      out
    )
  })
  tbl <- dplyr::bind_rows(rows) %>%
    dplyr::mutate(Term = ifelse(Term == "Intercept", Term, vapply(Term, friendly, character(1))))

  write_csv(tbl, file.path(OUT_DIR, paste0("coefficients_", resp_tag(resp_var), "_", model_name, ".csv")))
}

# -------------------- Effect plots + composite figure --------------------
make_effect_plot <- function(effect_spec) {
  parts <- strsplit(effect_spec, "->", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(ggplot() + theme_void())
  predictor <- parts[1]
  response <- parts[2]
  # response in effect_spec is variable name (e.g., temp, oxy, chl_log1p)
  if (!response %in% names(fit_piece)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste0("Missing fit for ", response), size = 5) + theme_void())
  }
  fit <- fit_piece[[response]]
  ce <- try(conditional_effects(fit, effects = predictor), silent = TRUE)
  if (inherits(ce, "try-error")) return(ggplot() + theme_void())
  pls <- try(plot(ce, plot = FALSE), silent = TRUE)
  if (inherits(pls, "try-error") || length(pls) == 0) return(ggplot() + theme_void())
  gp <- pls[[1]]

  # Add observed points from the equation's data (mf used in brm call)
  obs <- try(as.data.frame(model.frame(fit)), silent = TRUE)
  if (!inherits(obs, "try-error") && all(c(predictor, response) %in% names(obs))) {
    obs <- obs %>%
      dplyr::select(dplyr::all_of(c(predictor, response))) %>%
      dplyr::rename(x = dplyr::all_of(predictor), y = dplyr::all_of(response)) %>%
      dplyr::filter(is.finite(.data$x), is.finite(.data$y))
    gp <- gp + geom_point(data = obs, aes(x = .data$x, y = .data$y),
                          inherit.aes = FALSE, alpha = 0.35, size = 1.6, color = "black")
  }

  gp +
    labs(title = NULL, x = friendly(predictor), y = friendly(response)) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_blank()
    )
}

panel_label <- function(p, lbl) {
  p + annotate("text", x = -Inf, y = Inf, label = lbl,
               hjust = -0.15, vjust = 1.25, size = 7, fontface = "bold")
}

EFFECT_LIST <- if (length(EFFECT_PANELS)) EFFECT_PANELS else character()
if (length(EFFECT_LIST) < 2 && length(PATHS) > 0) {
  # derive
  term_of <- function(node) {
    switch(node,
      "chllog1p" = "chl_log1p",
      "callog1p" = "cal_log1p",
      "cyclog1p" = "cyc_log1p",
      node
    )
  }
  derived <- character()
  for (pth in PATHS) {
    if (length(pth) < 2) next
    for (k in seq_len(length(pth) - 1)) {
      derived <- c(derived, paste0(term_of(pth[k]), "->", term_of(pth[k + 1])))
    }
  }
  EFFECT_LIST <- unique(derived)
}

## -------------------- Path effects + standardized diagram --------------------
## Run the shared script *before* building the composite so the composite embeds
## the latest standardized `path_diagram.png` (and raw+std edge CSVs).
try({
  eff_script <- "sem_runner/path_effects_and_diagram_piecewise.R"
  if (file.exists(eff_script)) {
    Sys.setenv(SEM_OUT_DIR = OUT_DIR)
    Sys.setenv(SEM_SPEC_PATH = normalizePath(spec_path))
    Sys.unsetenv("SEM_FIT_RDS") # force piecewise backend when deep fit is absent
    message("Sourcing path effects script: ", eff_script)
    sys.source(eff_script, envir = environment())
  }
}, silent = FALSE)

# Build composite figure: path diagram rendered as raster in ggplot
png_path <- file.path(OUT_DIR, "path_diagram.png")
if (!file.exists(png_path)) stop("Expected path diagram not found: ", png_path)
imgA <- png::readPNG(png_path)
p_path <- ggplot() +
  annotation_custom(grid::rasterGrob(imgA, interpolate = TRUE), xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_fixed(clip = "off") +
  theme_void()

pA <- panel_label(p_path, "A")
p_effects <- lapply(EFFECT_LIST, make_effect_plot)
if (length(p_effects) > 0) {
  letters <- LETTERS[2:(1 + length(p_effects))]
  for (i in seq_along(p_effects)) p_effects[[i]] <- panel_label(p_effects[[i]], letters[i])
}

rows <- ceiling(length(p_effects) / 2)
effects_grid <- if (length(p_effects)) wrap_plots(p_effects, ncol = 2) else NULL
final_plot <- if (is.null(effects_grid)) (pA & theme(plot.background = element_blank())) else
  ((pA + theme(plot.margin = margin(0,0,6,0))) / (effects_grid + theme(plot.margin = margin(0,0,0,0)))) &
  theme(plot.background = element_blank(), plot.margin = margin(2,2,2,2))

base_h <- 6
per_row <- 3
plot_h <- if (length(p_effects) == 0) base_h else base_h + rows * per_row

ggsave(file.path(OUT_DIR, "composite_sem_AA_BC.tiff"),
       final_plot, width = 14, height = plot_h, dpi = 600, bg = "white", compression = "lzw")
ggsave(file.path(OUT_DIR, "composite_sem_AA_BC.png"),
       final_plot, width = 14, height = plot_h, dpi = 175, bg = "white")

# -------------------- README (reuse existing writer; add piecewise table) --------------------
final_lakes <- suppressWarnings(max(coverage_tbl$n_lakes, na.rm = TRUE))
if (!is.finite(final_lakes)) final_lakes <- NA_integer_
saveRDS(list(
  final_lakes = final_lakes,
  total_lakes_raw = total_lakes_raw,
  deep_single_candidates = deep_single_candidates,
  piecewise_n_by_equation = coverage_tbl,
  piecewise_missingness_summary = missingness_summary_tbl,
  piecewise_missingness_files = missingness_files
),
        file.path(OUT_DIR, "readme_context.rds"))

try({
  readme_script <- "sem_runner/write_run_readme.R"
  if (file.exists(readme_script)) {
    Sys.setenv(SEM_OUT_DIR = OUT_DIR)
    Sys.setenv(SEM_MODEL_NAME = model_name)
    Sys.setenv(SEM_SPEC_PATH = normalizePath(spec_path))
    # no SEM_FIT_RDS for piecewise
    message("Sourcing run README writer: ", readme_script)
    sys.source(readme_script, envir = environment())
  }
}, silent = TRUE)

message("=== SEM RUNNER (PIECEWISE BRMS; PER-EQUATION DATA): DONE ===")

