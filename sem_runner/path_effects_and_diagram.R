#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(posterior)
})

# Allow runner to pass current fit/output via env vars; fallback to defaults
OUT_DIR  <- Sys.getenv("SEM_OUT_DIR", unset = "sem_runner/results_deep")
FIT_ENV  <- Sys.getenv("SEM_FIT_RDS", unset = "")
PLOT_DIR <- OUT_DIR

deep_fit_path <- if (nzchar(FIT_ENV)) FIT_ENV else file.path(OUT_DIR, "sem_gaussian_log1p_rescor.rds")
piecewise_fit_paths <- list.files(OUT_DIR, pattern = "^fit_.*\\.rds$", full.names = TRUE)

use_deep <- file.exists(deep_fit_path)
use_piecewise <- length(piecewise_fit_paths) > 0

if (!use_deep && !use_piecewise) {
  stop("No SEM fit found. Expected either:\n",
       "- deep fit: ", deep_fit_path, "\n",
       "- piecewise fits: ", file.path(OUT_DIR, "fit_<response>.rds"))
}

fit <- NULL
ndr <- NULL
fits_pw <- list()
draws_pw <- list()

if (use_deep) {
  FIT_PATH <- deep_fit_path
  fit <- readRDS(FIT_PATH)
  ndr <- as_draws_df(fit)
} else {
  FIT_PATH <- ""
  for (fp in piecewise_fit_paths) {
    resp <- sub("^fit_(.+)\\.rds$", "\\1", basename(fp))
    fits_pw[[resp]] <- readRDS(fp)
    draws_pw[[resp]] <- posterior::as_draws_df(fits_pw[[resp]])
  }
}

# Optionally load spec to get explicit PATHS
SPEC_ENV <- Sys.getenv("SEM_SPEC_PATH", unset = "")
if (nzchar(SPEC_ENV) && file.exists(SPEC_ENV)) {
  try(sys.source(SPEC_ENV, envir = environment()), silent = TRUE)
}

# Friendly labels and mapping
friendly <- function(v) {
  v <- as.character(v)
  dplyr::case_when(
    v %in% c("chllog1p", "chl_log1p") ~ "chlorophyll",
    v %in% c("callog1p", "cal_log1p") ~ "calanoids",
    v %in% c("cyclog1p", "cyc_log1p") ~ "cyclopoids",
    v == "eutroph" ~ "eutrophication composite",
    v == "disease" ~ "copepod infection burden",
    v == "thermo" ~ "thermocline depth",
    v == "MaxDepth" ~ "max depth",
    v == "fdom" ~ "fDOM",
    v == "pH" ~ "pH",
    v == "oxy" ~ "oxygen",
    v == "temp" ~ "temperature",
    v == "defor" ~ "deforestation",
    TRUE ~ v
  )
}
resp_tag_of <- function(var) {
  switch(var,
    "chl_log1p" = "chllog1p",
    "cal_log1p" = "callog1p",
    "cyc_log1p" = "cyclog1p",
    var
  )
}

node_var_of <- function(node) {
  node <- as.character(node)
  dplyr::case_when(
    node == "chllog1p" ~ "chl_log1p",
    node == "callog1p" ~ "cal_log1p",
    node == "cyclog1p" ~ "cyc_log1p",
    TRUE ~ node
  )
}

# Detect present responses
candidate_resps <- c("chllog1p","callog1p","cyclog1p","pH","oxy","temp","eutroph","disease")
present_resps <- character()
if (use_deep) {
  present_resps <- candidate_resps[vapply(candidate_resps, function(r) {
    !inherits(try(fixef(fit, resp = r), silent = TRUE), "try-error")
  }, logical(1))]
} else {
  # piecewise: present responses are the available fit_<response>.rds files,
  # mapped to the response tags used in PATHS where relevant
  present_resps <- unique(c(names(fits_pw), vapply(names(fits_pw), resp_tag_of, character(1))))
}

# Build edges with term names for coefficient lookup (prefer explicit PATHS if provided)
edges <- list()
if (exists("PATHS", inherits = FALSE) && length(PATHS) > 0) {
  for (pth in PATHS) {
    if (length(pth) < 2) next
    for (k in seq_len(length(pth) - 1)) {
      from_node <- pth[k]
      to_node   <- pth[k + 1]
      if (use_deep) {
        if (!(to_node %in% present_resps)) next
        edges[[length(edges) + 1]] <- c(from = from_node, to = to_node, term = node_var_of(from_node))
      } else {
        from_v <- node_var_of(from_node)
        to_v   <- node_var_of(to_node)
        if (!to_v %in% names(draws_pw)) next
        edges[[length(edges) + 1]] <- c(from = from_v, to = to_v, term = from_v)
      }
    }
  }
} else {
  if (use_deep) {
    for (r in present_resps) {
      fx <- try(fixef(fit, resp = r), silent = TRUE)
      if (inherits(fx, "try-error")) next
      preds <- setdiff(rownames(fx), "Intercept")
      for (p in preds) edges[[length(edges) + 1]] <- c(from = resp_tag_of(p), to = r, term = p)
    }
  } else {
    for (to_v in names(draws_pw)) {
      fx <- try(fixef(fits_pw[[to_v]]), silent = TRUE)
      if (inherits(fx, "try-error")) next
      preds <- setdiff(rownames(fx), "Intercept")
      for (p in preds) edges[[length(edges) + 1]] <- c(from = p, to = to_v, term = p)
    }
  }
}
edges_df <- if (length(edges)) tibble::as_tibble(do.call(rbind, edges)) else tibble::tibble(from=character(), to=character(), term=character())

# Summary helpers
summ_df <- function(x) {
  x <- as.numeric(x)
  tibble(
    median = stats::median(x, na.rm = TRUE),
    q2.5 = as.numeric(stats::quantile(x, 0.025, na.rm = TRUE)),
    q97.5 = as.numeric(stats::quantile(x, 0.975, na.rm = TRUE)),
    n_draws = sum(is.finite(x))
  )
}
row_eff <- function(from, to, x) dplyr::bind_cols(tibble(path = paste(friendly(from), "->", friendly(to))), summ_df(x))

nonzero_95_flag <- function(q2.5, q97.5) {
  is.finite(q2.5) & is.finite(q97.5) & ((q2.5 > 0 & q97.5 > 0) | (q2.5 < 0 & q97.5 < 0))
}

# Formatting for edge labels (used in diagram + edge CSV)
fmt_edge <- function(x) {
  x <- as.numeric(x)
  if (!length(x)) return(NA_character_)
  m <- stats::median(x, na.rm = TRUE)
  if (!is.finite(m)) return(NA_character_)
  a <- abs(m)
  if (a == 0) return("0")
  # Use scientific notation for very small/large effects; otherwise show 3 decimals
  # to avoid different edges rounding to the same value (common after standardization).
  if (a < 1e-3 || a >= 100) return(formatC(m, format = "E", digits = 2))
  sprintf("%.3f", m)
}

# Standardization helpers ------------------------------------------------------
# We compute "standardized" effects by scaling an effect for X->Y by SD(X)/SD(Y).
# For indirect and total effects (which are already in units of Y per unit X),
# the same scaling yields the standardized effect.
source("sem_runner/scripts/sd_by_var_utils.R", local = TRUE)
sd_map_tbl <- sd_by_var_read(OUT_DIR)
if (is.null(sd_map_tbl)) {
  # Fallbacks when sd_by_var.csv is absent (e.g., deep runs or standalone calls)
  if (use_deep) {
    sd_map_tbl <- sd_by_var_compute(tryCatch(fit$data, error = function(e) NULL))
  } else {
    frames <- lapply(names(fits_pw), function(r) tryCatch(as.data.frame(model.frame(fits_pw[[r]])), error = function(e) NULL))
    frames <- frames[!vapply(frames, is.null, logical(1))]
    sd_map_tbl <- sd_by_var_compute_from_frames(frames)
  }
}
sd_lookup <- sd_by_var_lookup_factory(sd_map_tbl, resp_tag_of = resp_tag_of)

standardize_draws <- function(draws, from_var, to_var) {
  sx <- sd_lookup(from_var)
  sy <- sd_lookup(to_var)
  if (!is.finite(sx) || !is.finite(sy) || sy == 0) return(rep(NA_real_, length(draws)))
  as.numeric(draws) * (sx / sy)
}

# Canonicalize edges for path enumeration
edges_canon <- edges_df %>%
  dplyr::transmute(
    from = node_var_of(.data$from),
    to = node_var_of(.data$to),
    term = node_var_of(.data$term)
  ) %>%
  dplyr::distinct()

edge_key <- function(from, to) paste0(from, "->", to)
edge_present <- setNames(rep(TRUE, nrow(edges_canon)), edge_key(edges_canon$from, edges_canon$to))

get_edge_draws <- function(from_var, to_var) {
  if (use_deep) {
    to_tag <- resp_tag_of(to_var)
    nm <- paste0("b_", to_tag, "_", from_var)
    if (is.null(ndr) || !nm %in% names(ndr)) return(NULL)
    as.numeric(ndr[[nm]])
  } else {
    if (!to_var %in% names(draws_pw)) return(NULL)
    nm <- paste0("b_", from_var)
    if (!nm %in% names(draws_pw[[to_var]])) return(NULL)
    as.numeric(draws_pw[[to_var]][[nm]])
  }
}

trim_to <- function(x, n) {
  if (is.null(x)) return(NULL)
  x <- as.numeric(x)
  if (length(x) < n) return(x)
  x[seq_len(n)]
}

product_of_edges <- function(path_nodes) {
  # path_nodes is a character vector like c("defor","temp","oxy")
  if (length(path_nodes) < 2) return(NULL)
  draws_list <- list()
  for (k in seq_len(length(path_nodes) - 1)) {
    a <- path_nodes[k]
    b <- path_nodes[k + 1]
    dr <- get_edge_draws(a, b)
    if (is.null(dr)) return(NULL)
    draws_list[[length(draws_list) + 1]] <- dr
  }
  n <- min(vapply(draws_list, length, integer(1)))
  if (!is.finite(n) || n < 1) return(NULL)
  out <- rep(1, n)
  for (dr in draws_list) out <- out * trim_to(dr, n)
  out
}

all_paths <- function(from, to, adj, max_len = 10) {
  # DFS over a DAG-sized graph; returns list of node vectors
  res <- list()
  stack <- list(list(path = c(from)))
  while (length(stack)) {
    cur <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    p <- cur$path
    last <- p[[length(p)]]
    if (identical(last, to)) {
      res[[length(res) + 1]] <- p
      next
    }
    if (length(p) >= max_len) next
    nxt <- adj[[last]]
    if (is.null(nxt) || !length(nxt)) next
    for (v in nxt) {
      if (v %in% p) next
      stack[[length(stack) + 1]] <- list(path = c(p, v))
    }
  }
  res
}

nodes_all <- sort(unique(c(edges_canon$from, edges_canon$to)))
adj <- split(edges_canon$to, edges_canon$from)

# Direct effects (existing behavior, but keep)
direct_rows <- list()
for (i in seq_len(nrow(edges_df))) {
  to <- edges_df$to[i]; term <- edges_df$term[i]; from <- edges_df$from[i]
  if (use_deep) {
    nm <- paste0("b_", to, "_", term)
    if (!nm %in% names(ndr)) next
    direct_rows[[length(direct_rows) + 1]] <- row_eff(from, to, as.numeric(ndr[[nm]]))
  } else {
    if (!to %in% names(draws_pw)) next
    nm <- paste0("b_", term)
    if (!nm %in% names(draws_pw[[to]])) next
    direct_rows[[length(direct_rows) + 1]] <- row_eff(from, to, as.numeric(draws_pw[[to]][[nm]]))
  }
}
effects_tbl <- if (length(direct_rows)) dplyr::bind_rows(direct_rows) else tibble::tibble(path=character(), median=numeric(), q2.5=numeric(), q97.5=numeric())
effects_tbl <- effects_tbl %>%
  dplyr::mutate(nonzero_95 = nonzero_95_flag(.data$q2.5, .data$q97.5))

# Add standardized direct effects (using the parsed from/to from the path string)
if (nrow(effects_tbl)) {
  # parse: "<from label>-><to label>" is NOT safe; re-use edges_canon to compute std draws
  std_rows <- list()
  for (i in seq_len(nrow(edges_canon))) {
    a <- edges_canon$from[i]; b <- edges_canon$to[i]
    dr <- get_edge_draws(a, b)
    if (is.null(dr)) next
    drs <- standardize_draws(dr, a, b)
    std_rows[[length(std_rows) + 1]] <- tibble::tibble(
      path = paste(friendly(a), "->", friendly(b))
    ) %>% dplyr::bind_cols(
      summ_df(drs) %>%
        dplyr::rename(median_std = .data$median, q2.5_std = .data$q2.5, q97.5_std = .data$q97.5, n_draws_std = .data$n_draws)
    )
  }
  std_tbl <- if (length(std_rows)) dplyr::bind_rows(std_rows) else tibble::tibble(path=character())
  effects_tbl <- effects_tbl %>%
    dplyr::left_join(std_tbl, by = "path") %>%
    dplyr::mutate(nonzero_95_std = nonzero_95_flag(.data$q2.5_std, .data$q97.5_std))
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

out_csv <- file.path(OUT_DIR, "path_effects.csv")
readr::write_csv(effects_tbl, out_csv)
message("Wrote ", out_csv)

# -------------------- Indirect + total effects (draw-by-draw products) --------------------
ind_rows <- list()
tot_rows <- list()

max_len <- length(nodes_all) + 1
for (a in nodes_all) {
  for (b in nodes_all) {
    if (identical(a, b)) next
    # only consider reachable pairs
    ps <- all_paths(a, b, adj, max_len = max_len)
    if (!length(ps)) next

    # Direct draws if direct edge exists
    direct_dr <- if (edge_key(a, b) %in% names(edge_present)) get_edge_draws(a, b) else NULL

    # Indirect paths: length >= 3 nodes
    indirect_ps <- ps[vapply(ps, length, integer(1)) >= 3]
    indirect_draws <- list()
    if (length(indirect_ps)) {
      for (j in seq_along(indirect_ps)) {
        p <- indirect_ps[[j]]
        dr <- product_of_edges(p)
        if (is.null(dr)) next
        indirect_draws[[length(indirect_draws) + 1]] <- dr
        via <- if (length(p) > 2) paste(p[2:(length(p) - 1)], collapse = " -> ") else ""
        ind_rows[[length(ind_rows) + 1]] <- tibble::tibble(
          from = a,
          to = b,
          via = via,
          path = paste(p, collapse = " -> ")
        ) %>% dplyr::bind_cols(summ_df(dr))
      }
    }

    # Total = direct + sum(indirect paths) draw-by-draw
    comp <- c(list(direct_dr), indirect_draws)
    comp <- comp[!vapply(comp, is.null, logical(1))]
    if (!length(comp)) next
    n <- min(vapply(comp, length, integer(1)))
    if (!is.finite(n) || n < 1) next
    total_dr <- rep(0, n)
    if (!is.null(direct_dr)) total_dr <- total_dr + trim_to(direct_dr, n)
    if (length(indirect_draws)) {
      for (dr in indirect_draws) total_dr <- total_dr + trim_to(dr, n)
    }

    tot_rows[[length(tot_rows) + 1]] <- tibble::tibble(
      from = a,
      to = b,
      direct_included = !is.null(direct_dr),
      n_indirect_paths = length(indirect_draws)
    ) %>% dplyr::bind_cols(summ_df(total_dr))
  }
}

ind_tbl <- if (length(ind_rows)) dplyr::bind_rows(ind_rows) else tibble::tibble(from=character(),to=character(),via=character(),path=character(),median=numeric(),q2.5=numeric(),q97.5=numeric(),n_draws=integer())
tot_tbl <- if (length(tot_rows)) dplyr::bind_rows(tot_rows) else tibble::tibble(from=character(),to=character(),direct_included=logical(),n_indirect_paths=integer(),median=numeric(),q2.5=numeric(),q97.5=numeric(),n_draws=integer())

# Add friendly labels for readability
ind_tbl <- ind_tbl %>%
  dplyr::mutate(nonzero_95 = nonzero_95_flag(.data$q2.5, .data$q97.5)) %>%
  dplyr::mutate(from_label = friendly(.data$from), to_label = friendly(.data$to))
tot_tbl <- tot_tbl %>%
  dplyr::mutate(nonzero_95 = nonzero_95_flag(.data$q2.5, .data$q97.5)) %>%
  dplyr::mutate(from_label = friendly(.data$from), to_label = friendly(.data$to))

# Standardize indirect/total effects
if (nrow(ind_tbl)) {
  ind_tbl <- ind_tbl %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      median_std = median * (sd_lookup(from) / sd_lookup(to)),
      q2.5_std = q2.5 * (sd_lookup(from) / sd_lookup(to)),
      q97.5_std = q97.5 * (sd_lookup(from) / sd_lookup(to)),
      nonzero_95_std = nonzero_95_flag(q2.5_std, q97.5_std)
    ) %>%
    dplyr::ungroup()
}
if (nrow(tot_tbl)) {
  tot_tbl <- tot_tbl %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      median_std = median * (sd_lookup(from) / sd_lookup(to)),
      q2.5_std = q2.5 * (sd_lookup(from) / sd_lookup(to)),
      q97.5_std = q97.5 * (sd_lookup(from) / sd_lookup(to)),
      nonzero_95_std = nonzero_95_flag(q2.5_std, q97.5_std)
    ) %>%
    dplyr::ungroup()
}

out_ind <- file.path(OUT_DIR, "path_effects_indirect.csv")
out_tot <- file.path(OUT_DIR, "path_effects_total.csv")
readr::write_csv(ind_tbl, out_ind)
readr::write_csv(tot_tbl, out_tot)
message("Wrote ", out_ind)
message("Wrote ", out_tot)

# Edge table backing the diagram labels (raw + standardized)
edge_rows <- list()
for (i in seq_len(nrow(edges_canon))) {
  a <- edges_canon$from[i]
  b <- edges_canon$to[i]
  dr_raw <- get_edge_draws(a, b)
  if (is.null(dr_raw)) next
  dr_std <- standardize_draws(dr_raw, a, b)
  edge_rows[[length(edge_rows) + 1]] <- tibble::tibble(
    from = a,
    to = b,
    from_label = friendly(a),
    to_label = friendly(b),
    median_raw = stats::median(dr_raw, na.rm = TRUE),
    q2.5_raw = as.numeric(stats::quantile(dr_raw, 0.025, na.rm = TRUE)),
    q97.5_raw = as.numeric(stats::quantile(dr_raw, 0.975, na.rm = TRUE)),
    nonzero_95_raw = nonzero_95_flag(q2.5_raw, q97.5_raw),
    median_std = stats::median(dr_std, na.rm = TRUE),
    q2.5_std = as.numeric(stats::quantile(dr_std, 0.025, na.rm = TRUE)),
    q97.5_std = as.numeric(stats::quantile(dr_std, 0.975, na.rm = TRUE)),
    nonzero_95_std = nonzero_95_flag(q2.5_std, q97.5_std),
    label_std = fmt_edge(dr_std)
  )
}
edge_tbl <- if (length(edge_rows)) dplyr::bind_rows(edge_rows) else tibble::tibble()
out_edge <- file.path(OUT_DIR, "path_diagram_edges.csv")
readr::write_csv(edge_tbl, out_edge)
message("Wrote ", out_edge)

# Path diagram (fail-fast if dependencies are missing; fallbacks are disabled)
require_pkg <- function(pkg, purpose) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required %s. Install it and re-run the SEM runner.", pkg, purpose), call. = FALSE)
  }
}

require_pkg("DiagrammeR", "to build the SEM path diagram")
require_pkg("DiagrammeRsvg", "to export the diagram as SVG")
has_magick <- requireNamespace("magick", quietly = TRUE)
has_rsvg   <- requireNamespace("rsvg", quietly = TRUE)
has_png    <- requireNamespace("png", quietly = TRUE)
if (!has_magick && !(has_rsvg && has_png)) {
  stop("Either the 'magick' package or both 'rsvg' and 'png' are required to render the SEM path diagram. Install the missing packages and re-run the SEM runner.", call. = FALSE)
}

library(DiagrammeR)
LABEL_MODE <- "std"
# nodes from edges or present responses if none
nod_ids <- unique(c(edges_df$from, edges_df$to))
if (!length(nod_ids)) nod_ids <- present_resps
nodes <- setNames(vapply(nod_ids, friendly, character(1)), nod_ids)
# edges with labels (direct medians)
edge_lines <- c()
for (i in seq_len(nrow(edges_df))) {
  to <- edges_df$to[i]; term <- edges_df$term[i]; from <- edges_df$from[i]
  if (use_deep) {
    nm <- paste0("b_", to, "_", term)
    if (!nm %in% names(ndr)) next
    from_v <- node_var_of(from)
    to_v   <- node_var_of(to)
    lab_draws <- standardize_draws(ndr[[nm]], from_v, to_v)
    edge_lines <- c(edge_lines, sprintf("%s -> %s [label=\"%s\"];", from_v, to_v, fmt_edge(lab_draws)))
  } else {
    nm <- paste0("b_", term)
    if (!to %in% names(draws_pw)) next
    if (!nm %in% names(draws_pw[[to]])) next
    from_v <- node_var_of(from)
    to_v   <- node_var_of(to)
    lab_draws <- standardize_draws(draws_pw[[to]][[nm]], from_v, to_v)
    edge_lines <- c(edge_lines, sprintf("%s -> %s [label=\"%s\"];", from_v, to_v, fmt_edge(lab_draws)))
  }
}
node_decl <- if (length(nodes)) paste(sprintf("%s [label=\"%s\"];", names(nodes), nodes), collapse = " ") else ""
edges_decl <- paste(edge_lines, collapse = "\n  ")
dot <- sprintf('digraph sem {\n  graph [rankdir=LR, margin=\"0\", pad=\"0\", bgcolor=\"#FFFFFF\"]\n  node [shape=box, style=filled, fillcolor=\"#FFFFFF\", color=\"black\"]\n  %s\n  %s\n}', node_decl, edges_decl)

png_path <- file.path(PLOT_DIR, "path_diagram.png")
svg_txt <- DiagrammeRsvg::export_svg(DiagrammeR::grViz(dot))
if (has_magick) {
  img <- magick::image_read_svg(svg_txt, width = 1200, height = 600)
  # Ensure a white background even if SVG retains transparency
  img <- magick::image_background(img, "white")
  magick::image_write(img, path = png_path, format = "png")
} else {
  rsvg::rsvg_png(charToRaw(svg_txt), file = png_path, width = 1200, height = 600)
}
message("Wrote ", png_path)
