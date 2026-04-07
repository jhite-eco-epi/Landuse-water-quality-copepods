#!/usr/bin/env Rscript
#
# Piecewise-only path effects + diagram builder.
# Assumes the run directory contains:
# - fit_<response>.rds files (one brmsfit per equation)
# - sd_by_var.csv (global SDs computed from deep/single lake medians)
# - SEM_SPEC_PATH pointing to a spec with PATHS (preferred) and/or EFFECT_PANELS
#
# Writes:
# - path_effects.csv (direct)
# - path_effects_indirect.csv (per indirect path)
# - path_effects_total.csv (total effects)
# - path_diagram_edges.csv (raw + standardized edge summaries)
# - path_diagram.png (standardized edge labels)

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(posterior)
})

OUT_DIR  <- Sys.getenv("SEM_OUT_DIR", unset = "sem_runner/results_piecewise")
SPEC_ENV <- Sys.getenv("SEM_SPEC_PATH", unset = "")

if (!dir.exists(OUT_DIR)) stop("Output directory not found: ", OUT_DIR)

fit_paths <- list.files(OUT_DIR, pattern = "^fit_.*\\.rds$", full.names = TRUE)
if (!length(fit_paths)) stop("No piecewise fits found in OUT_DIR (expected fit_<response>.rds).")

fits <- list()
draws <- list()
for (fp in fit_paths) {
  resp <- sub("^fit_(.+)\\.rds$", "\\1", basename(fp))
  fits[[resp]] <- readRDS(fp)
  draws[[resp]] <- posterior::as_draws_df(fits[[resp]])
}

# Load spec for PATHS
if (nzchar(SPEC_ENV) && file.exists(SPEC_ENV)) {
  try(sys.source(SPEC_ENV, envir = environment()), silent = TRUE)
}
if (!exists("PATHS", inherits = FALSE) || length(PATHS) == 0) {
  stop("PATHS not found in spec (SEM_SPEC_PATH). Piecewise effects require PATHS to define the graph.")
}

source("sem_runner/scripts/sd_by_var_utils.R", local = TRUE)
sd_tbl <- sd_by_var_read(OUT_DIR)
if (is.null(sd_tbl) || !nrow(sd_tbl)) stop("Missing or empty sd_by_var.csv in OUT_DIR.")
sd_lookup <- sd_by_var_lookup_factory(sd_tbl)

standardize_draws <- function(draws, from_var, to_var) {
  sx <- sd_lookup(from_var)
  sy <- sd_lookup(to_var)
  if (!is.finite(sx) || !is.finite(sy) || sy == 0) return(rep(NA_real_, length(draws)))
  as.numeric(draws) * (sx / sy)
}

nonzero_95_flag <- function(q2.5, q97.5) {
  is.finite(q2.5) & is.finite(q97.5) & ((q2.5 > 0 & q97.5 > 0) | (q2.5 < 0 & q97.5 < 0))
}

summ_draws <- function(x) {
  x <- as.numeric(x)
  tibble(
    median = stats::median(x, na.rm = TRUE),
    q2.5 = as.numeric(stats::quantile(x, 0.025, na.rm = TRUE)),
    q97.5 = as.numeric(stats::quantile(x, 0.975, na.rm = TRUE)),
    n_draws = sum(is.finite(x))
  )
}

fmt_edge <- function(x) {
  x <- as.numeric(x)
  if (!length(x)) return(NA_character_)
  m <- stats::median(x, na.rm = TRUE)
  if (!is.finite(m)) return(NA_character_)
  a <- abs(m)
  if (a == 0) return("0")
  if (a < 1e-3 || a >= 100) return(formatC(m, format = "E", digits = 2))
  sprintf("%.3f", m)
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

friendly <- function(v) {
  v <- as.character(v)
  dplyr::case_when(
    v %in% c("chllog1p", "chl_log1p") ~ "chlorophyll",
    v %in% c("callog1p", "cal_log1p") ~ "calanoids",
    v %in% c("cyclog1p", "cyc_log1p") ~ "cyclopoids",
    v == "eutroph" ~ "eutrophication composite",
    v == "thermo" ~ "thermocline depth",
    v == "pH" ~ "pH",
    v == "oxy" ~ "oxygen",
    v == "temp" ~ "temperature",
    v == "defor" ~ "deforestation",
    TRUE ~ v
  )
}

# Build canonical edge list from PATHS
edges <- list()
for (pth in PATHS) {
  if (length(pth) < 2) next
  for (k in seq_len(length(pth) - 1)) {
    a <- node_var_of(pth[k])
    b <- node_var_of(pth[k + 1])
    edges[[length(edges) + 1]] <- c(from = a, to = b)
  }
}
edges_df <- if (length(edges)) tibble::as_tibble(do.call(rbind, edges)) else tibble(from = character(), to = character())
edges_df <- edges_df %>% dplyr::distinct()

get_edge_draws <- function(from_var, to_var) {
  if (!to_var %in% names(draws)) return(NULL)
  nm <- paste0("b_", from_var)
  if (!nm %in% names(draws[[to_var]])) return(NULL)
  as.numeric(draws[[to_var]][[nm]])
}

trim_to <- function(x, n) {
  x <- as.numeric(x)
  if (length(x) < n) return(x)
  x[seq_len(n)]
}

product_of_edges <- function(path_nodes) {
  if (length(path_nodes) < 2) return(NULL)
  ds <- list()
  for (k in seq_len(length(path_nodes) - 1)) {
    a <- path_nodes[k]; b <- path_nodes[k + 1]
    dr <- get_edge_draws(a, b)
    if (is.null(dr)) return(NULL)
    ds[[length(ds) + 1]] <- dr
  }
  n <- min(vapply(ds, length, integer(1)))
  if (!is.finite(n) || n < 1) return(NULL)
  out <- rep(1, n)
  for (dr in ds) out <- out * trim_to(dr, n)
  out
}

all_paths <- function(from, to, adj, max_len = 20) {
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

nodes_all <- sort(unique(c(edges_df$from, edges_df$to)))
adj <- split(edges_df$to, edges_df$from)
edge_present <- setNames(rep(TRUE, nrow(edges_df)), paste0(edges_df$from, "->", edges_df$to))

# Direct effects
direct_rows <- list()
for (i in seq_len(nrow(edges_df))) {
  a <- edges_df$from[i]; b <- edges_df$to[i]
  dr <- get_edge_draws(a, b)
  if (is.null(dr)) next
  drs <- standardize_draws(dr, a, b)
  direct_rows[[length(direct_rows) + 1]] <- tibble::tibble(
    path = paste(friendly(a), "->", friendly(b))
  ) %>% dplyr::bind_cols(summ_draws(dr)) %>%
    dplyr::mutate(nonzero_95 = nonzero_95_flag(.data$q2.5, .data$q97.5)) %>%
    dplyr::bind_cols(
      summ_draws(drs) %>%
        dplyr::rename(median_std = .data$median, q2.5_std = .data$q2.5, q97.5_std = .data$q97.5, n_draws_std = .data$n_draws)
    ) %>%
    dplyr::mutate(nonzero_95_std = nonzero_95_flag(.data$q2.5_std, .data$q97.5_std))
}
direct_tbl <- if (length(direct_rows)) dplyr::bind_rows(direct_rows) else tibble::tibble()
readr::write_csv(direct_tbl, file.path(OUT_DIR, "path_effects.csv"))

# Indirect + total effects
ind_rows <- list()
tot_rows <- list()
for (a in nodes_all) {
  for (b in nodes_all) {
    if (identical(a, b)) next
    ps <- all_paths(a, b, adj, max_len = length(nodes_all) + 1)
    if (!length(ps)) next

    direct_dr <- if (paste0(a, "->", b) %in% names(edge_present)) get_edge_draws(a, b) else NULL
    indirect_ps <- ps[vapply(ps, length, integer(1)) >= 3]
    indirect_draws <- list()
    if (length(indirect_ps)) {
      for (p in indirect_ps) {
        dr <- product_of_edges(p)
        if (is.null(dr)) next
        indirect_draws[[length(indirect_draws) + 1]] <- dr
        via <- paste(p[2:(length(p) - 1)], collapse = " -> ")
        ind_rows[[length(ind_rows) + 1]] <- tibble::tibble(
          from = a,
          to = b,
          via = via,
          path = paste(p, collapse = " -> ")
        ) %>% dplyr::bind_cols(summ_draws(dr)) %>%
          dplyr::mutate(nonzero_95 = nonzero_95_flag(.data$q2.5, .data$q97.5)) %>%
          dplyr::mutate(
            median_std = median * (sd_lookup(a) / sd_lookup(b)),
            q2.5_std = q2.5 * (sd_lookup(a) / sd_lookup(b)),
            q97.5_std = q97.5 * (sd_lookup(a) / sd_lookup(b)),
            nonzero_95_std = nonzero_95_flag(.data$q2.5_std, .data$q97.5_std)
          ) %>%
          dplyr::mutate(from_label = friendly(a), to_label = friendly(b))
      }
    }

    comp <- c(list(direct_dr), indirect_draws)
    comp <- comp[!vapply(comp, is.null, logical(1))]
    if (!length(comp)) next
    n <- min(vapply(comp, length, integer(1)))
    if (!is.finite(n) || n < 1) next
    total_dr <- rep(0, n)
    if (!is.null(direct_dr)) total_dr <- total_dr + trim_to(direct_dr, n)
    for (dr in indirect_draws) total_dr <- total_dr + trim_to(dr, n)

    s <- summ_draws(total_dr)
    tot_rows[[length(tot_rows) + 1]] <- tibble::tibble(
      from = a, to = b,
      direct_included = !is.null(direct_dr),
      n_indirect_paths = length(indirect_draws)
    ) %>%
      dplyr::bind_cols(s) %>%
      dplyr::mutate(nonzero_95 = nonzero_95_flag(.data$q2.5, .data$q97.5)) %>%
      dplyr::mutate(
        median_std = median * (sd_lookup(a) / sd_lookup(b)),
        q2.5_std = q2.5 * (sd_lookup(a) / sd_lookup(b)),
        q97.5_std = q97.5 * (sd_lookup(a) / sd_lookup(b)),
        nonzero_95_std = nonzero_95_flag(.data$q2.5_std, .data$q97.5_std),
        from_label = friendly(a),
        to_label = friendly(b)
      )
  }
}

ind_tbl <- if (length(ind_rows)) dplyr::bind_rows(ind_rows) else tibble::tibble()
tot_tbl <- if (length(tot_rows)) dplyr::bind_rows(tot_rows) else tibble::tibble()
readr::write_csv(ind_tbl, file.path(OUT_DIR, "path_effects_indirect.csv"))
readr::write_csv(tot_tbl, file.path(OUT_DIR, "path_effects_total.csv"))

# Edge table backing diagram labels
edge_rows <- list()
for (i in seq_len(nrow(edges_df))) {
  a <- edges_df$from[i]; b <- edges_df$to[i]
  dr_raw <- get_edge_draws(a, b)
  if (is.null(dr_raw)) next
  dr_std <- standardize_draws(dr_raw, a, b)
  edge_rows[[length(edge_rows) + 1]] <- tibble::tibble(
    from = a, to = b,
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
readr::write_csv(edge_tbl, file.path(OUT_DIR, "path_diagram_edges.csv"))

# Diagram PNG (standardized labels only)
require_pkg <- function(pkg, purpose) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required %s.", pkg, purpose), call. = FALSE)
  }
}
require_pkg("DiagrammeR", "to build the SEM path diagram")
require_pkg("DiagrammeRsvg", "to export the diagram as SVG")
has_magick <- requireNamespace("magick", quietly = TRUE)
has_rsvg   <- requireNamespace("rsvg", quietly = TRUE)
has_png    <- requireNamespace("png", quietly = TRUE)
if (!has_magick && !(has_rsvg && has_png)) {
  stop("Either the 'magick' package or both 'rsvg' and 'png' are required to render the SEM path diagram.", call. = FALSE)
}

node_ids <- unique(c(edges_df$from, edges_df$to))
nodes <- setNames(vapply(node_ids, friendly, character(1)), node_ids)
node_decl <- paste(sprintf("%s [label=\"%s\"];", names(nodes), nodes), collapse = " ")

edge_lines <- c()
for (i in seq_len(nrow(edge_tbl))) {
  edge_lines <- c(edge_lines, sprintf("%s -> %s [label=\"%s\"];", edge_tbl$from[i], edge_tbl$to[i], edge_tbl$label_std[i]))
}
dot <- sprintf('digraph sem {\n  graph [rankdir=LR, margin=\"0\", pad=\"0\", bgcolor=\"#FFFFFF\"]\n  node [shape=box, style=filled, fillcolor=\"#FFFFFF\", color=\"black\"]\n  %s\n  %s\n}', node_decl, paste(edge_lines, collapse = "\n  "))

svg_txt <- DiagrammeRsvg::export_svg(DiagrammeR::grViz(dot))
png_path <- file.path(OUT_DIR, "path_diagram.png")
if (has_magick) {
  img <- magick::image_read_svg(svg_txt, width = 1200, height = 600)
  img <- magick::image_background(img, "white")
  magick::image_write(img, path = png_path, format = "png")
} else {
  rsvg::rsvg_png(charToRaw(svg_txt), file = png_path, width = 1200, height = 600)
}

message("Wrote path effects and diagram to ", OUT_DIR)

