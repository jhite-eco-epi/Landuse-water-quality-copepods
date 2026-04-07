#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

DATA_RDS_DEFAULT <- "data/main_zooplankton_data.rds"

norm_id <- function(x) tolower(trimws(as.character(x)))
safe_num <- function(x) suppressWarnings(as.numeric(x))

load_main <- function(data_rds = DATA_RDS_DEFAULT) {
  if (!file.exists(data_rds)) stop("Missing main RDS: ", data_rds)
  readRDS(data_rds)
}

source_postprocessing <- function() {
  source("data_pipeline/postprocessing.R")
}

filter_deep_single_as_deep <- function(df) {
  if (!"zoop_sample_type" %in% names(df)) stop("Missing zoop_sample_type in main.")
  df %>%
    mutate(zoop_sample_type = tolower(trimws(as.character(zoop_sample_type)))) %>%
    filter(zoop_sample_type %in% c("deep", "single")) %>%
    mutate(zoop_sample_type = "deep")
}

lake_medians <- function(df) {
  source_postprocessing()
  if (!"LakeID" %in% names(df)) stop("Missing LakeID in main.")
  df %>%
    mutate(LakeID = norm_id(LakeID)) %>%
    group_by(LakeID) %>%
    summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)), .groups = "drop") %>%
    exclude_lakes(.warn = FALSE)
}

pretty_taxon <- function(x) {
  x <- sub("_density_per_L$", "", x)
  x <- gsub("\\.+", " ", x)
  x <- gsub("_+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

order_correlation_matrix <- function(cor_mat, order = NULL, cluster_by_abs = FALSE) {
  if (!is.null(order)) {
    keep <- order[order %in% colnames(cor_mat)]
    cor_mat[keep, keep, drop = FALSE]
  } else if (isTRUE(cluster_by_abs)) {
    dist_mat <- as.dist(1 - abs(cor_mat))
    ord <- hclust(dist_mat)$order
    cor_mat[ord, ord, drop = FALSE]
  } else {
    cor_mat
  }
}

correlation_triangle_df <- function(cor_mat, labels_map = NULL, triangle = c("lower", "upper"), digits = 2) {
  triangle <- match.arg(triangle)
  var_order <- colnames(cor_mat)
  pretty_labels <- if (is.null(labels_map)) {
    setNames(var_order, var_order)
  } else {
    labels_map[var_order]
  }

  as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE) %>%
    rename(var1 = Var1, var2 = Var2, corr = Freq) %>%
    mutate(
      var1_raw = var1,
      var2_raw = var2,
      i = match(var1_raw, var_order),
      j = match(var2_raw, var_order)
    ) %>%
    filter(if (triangle == "lower") j >= i else j <= i) %>%
    mutate(
      var1 = factor(var1_raw, levels = var_order, labels = unname(pretty_labels[var_order])),
      var2 = factor(var2_raw, levels = rev(var_order), labels = unname(pretty_labels[rev(var_order)])),
      label = sprintf(paste0("%.", digits, "f"), corr)
    )
}

plot_correlation_heatmap <- function(
  cor_df,
  fill_name = "Correlation",
  low = "darkcyan",
  mid = "white",
  high = "#8b008b",
  midpoint = 0,
  limits = c(-1, 1),
  bold_threshold = 0.7,
  text_size = 3.4
) {
  ggplot(cor_df, aes(x = var1, y = var2, fill = corr)) +
    geom_tile(color = "grey92") +
    geom_text(
      aes(label = label, fontface = ifelse(abs(corr) >= bold_threshold, "bold", "plain")),
      size = text_size
    ) +
    scale_fill_gradient2(
      low = low,
      mid = mid,
      high = high,
      midpoint = midpoint,
      limits = limits,
      name = fill_name
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_blank()
    )
}

correlation_table <- function(cor_df, value_name = "correlation") {
  out_tbl <- cor_df %>%
    transmute(
      variable1 = var1_raw,
      variable2 = var2_raw,
      label1 = as.character(var1),
      label2 = as.character(var2),
      !!value_name := corr
    )
  out_tbl[order(-abs(out_tbl[[value_name]])), , drop = FALSE]
}

save_correlation_outputs <- function(
  plot,
  table,
  out_dir,
  stats_dir = out_dir,
  filename_base,
  table_filename,
  width = 12,
  height = 11
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

  png_path <- file.path(out_dir, paste0(filename_base, ".png"))
  tiff_path <- file.path(out_dir, paste0(filename_base, ".tiff"))
  table_path <- file.path(stats_dir, table_filename)

  ggsave(png_path, plot, width = width, height = height, dpi = 200, bg = "white")
  ggsave(tiff_path, plot, width = width, height = height, dpi = 600, bg = "white", compression = "lzw")
  readr::write_csv(table, table_path)

  invisible(list(png = png_path, tiff = tiff_path, table = table_path))
}

write_correlations_results_readme <- function(
  out_path = "correlations/outputs/README.md",
  biotic_png = "biotic/biotic_correlations.png",
  biotic_table = "biotic/biotic_correlation_table.csv",
  abiotic_png = "abiotic/abiotic_correlations.png",
  abiotic_table = "abiotic/abiotic_correlation_table.csv"
) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  lines <- c(
    "# Correlation Results",
    "",
    "This report shows the current correlation figures generated from the shared analysis dataset.",
    "",
    "## Biotic Correlations",
    "",
    "Lake-level Spearman correlations across zooplankton taxa.",
    "",
    paste0("![Biotic correlations](", biotic_png, ")"),
    "",
    paste0("- Table: `", biotic_table, "`"),
    "",
    "## Abiotic Correlations",
    "",
    "Lake-level Spearman correlations across integrated environmental predictors and thermocline depth.",
    "",
    paste0("![Abiotic correlations](", abiotic_png, ")"),
    "",
    paste0("- Table: `", abiotic_table, "`")
  )

  writeLines(lines, con = out_path)
  invisible(out_path)
}

prepare_biotic_correlation_result <- function(
  data_rds = DATA_RDS_DEFAULT,
  min_obs = 6
) {
  df <- load_main(data_rds) %>%
    filter_deep_single_as_deep()

  d <- lake_medians(df)

  taxa_cols <- names(d)[grepl("_density_per_L$", names(d))]
  if (!length(taxa_cols)) stop("No *_density_per_L columns found after aggregation.")

  taxa_pretty_map <- setNames(vapply(taxa_cols, pretty_taxon, character(1)), taxa_cols)

  taxa_df <- d %>%
    dplyr::select(LakeID, all_of(taxa_cols)) %>%
    dplyr::mutate(dplyr::across(all_of(taxa_cols), safe_num))

  keep <- vapply(taxa_cols, function(cc) {
    x <- taxa_df[[cc]]
    n_ok <- sum(is.finite(x))
    sd_ok <- suppressWarnings(sd(x, na.rm = TRUE))
    is.finite(n_ok) && n_ok >= min_obs && is.finite(sd_ok) && sd_ok > 0
  }, logical(1))

  kept_cols <- taxa_cols[keep]
  if (length(kept_cols) < 3) stop("Too few taxa columns after filtering (need >= 3).")

  mat <- as.matrix(taxa_df[, kept_cols, drop = FALSE])
  cor_mat <- suppressWarnings(cor(mat, use = "pairwise.complete.obs", method = "spearman"))
  cor_mat <- order_correlation_matrix(cor_mat, cluster_by_abs = TRUE)

  cor_df <- correlation_triangle_df(
    cor_mat = cor_mat,
    labels_map = taxa_pretty_map,
    triangle = "lower",
    digits = 2
  )

  plot <- plot_correlation_heatmap(
    cor_df = cor_df,
    fill_name = "Spearman r",
    low = "darkcyan",
    mid = "white",
    high = "#8b008b",
    bold_threshold = 0.7,
    text_size = 3.4
  )

  table <- correlation_table(cor_df, value_name = "spearman_r")

  list(plot = plot, cor_mat = cor_mat, table = table, plot_df = cor_df)
}

build_biotic_correlations <- function(
  out_dir = "correlations/outputs/biotic",
  stats_dir = out_dir,
  data_rds = DATA_RDS_DEFAULT,
  filename_base = "biotic_correlations"
) {
  result <- prepare_biotic_correlation_result(data_rds = data_rds)
  files <- save_correlation_outputs(
    plot = result$plot,
    table = result$table,
    out_dir = out_dir,
    stats_dir = stats_dir,
    filename_base = filename_base,
    table_filename = "biotic_correlation_table.csv",
    width = 12,
    height = 11
  )
  result$files <- files
  invisible(result)
}

prepare_abiotic_correlation_result <- function(
  data_rds = DATA_RDS_DEFAULT,
  min_depth_m = 2
) {
  df <- load_main(data_rds) %>%
    filter_deep_single_as_deep() %>%
    mutate(LakeID = norm_id(LakeID))

  lake_df <- df %>%
    group_by(LakeID) %>%
    summarise(
      integrated_chl = median(safe_num(integrated_chl), na.rm = TRUE),
      integrated_fDOM = median(safe_num(integrated_fDOM), na.rm = TRUE),
      integrated_pH = median(safe_num(integrated_pH), na.rm = TRUE),
      integrated_temp = median(safe_num(integrated_temp), na.rm = TRUE),
      integrated_DO_percent = median(safe_num(integrated_DO_percent), na.rm = TRUE),
      thermocline_depth_m = median(safe_num(thermocline_depth_m), na.rm = TRUE),
      MaxDepth = median(safe_num(MaxDepth), na.rm = TRUE),
      .groups = "drop"
    )

  source_postprocessing()
  lake_df <- exclude_lakes(lake_df, .warn = FALSE)
  lake_df <- lake_df %>%
    filter(!is.finite(MaxDepth) | MaxDepth >= min_depth_m)

  vars <- c(
    "integrated_chl",
    "integrated_fDOM",
    "integrated_pH",
    "integrated_temp",
    "integrated_DO_percent",
    "thermocline_depth_m"
  )

  missing_vars <- vars[!vars %in% names(lake_df)]
  if (length(missing_vars)) stop("Missing abiotic columns: ", paste(missing_vars, collapse = ", "))

  abiotic_df <- lake_df %>%
    dplyr::select(all_of(vars))

  labels_map <- c(
    integrated_chl = "Chlorophyll",
    integrated_fDOM = "fDOM",
    integrated_pH = "pH",
    integrated_temp = "Temperature",
    integrated_DO_percent = "DO",
    thermocline_depth_m = "Thermocline Depth"
  )

  cor_mat <- cor(abiotic_df, use = "pairwise.complete.obs", method = "spearman")
  cor_mat <- order_correlation_matrix(cor_mat, order = vars)

  cor_df <- correlation_triangle_df(
    cor_mat = cor_mat,
    labels_map = labels_map,
    triangle = "lower",
    digits = 2
  )

  plot <- plot_correlation_heatmap(
    cor_df = cor_df,
    fill_name = "Spearman r",
    low = "darkcyan",
    mid = "white",
    high = "#8b008b",
    bold_threshold = 0.55,
    text_size = 4
  )

  table <- correlation_table(cor_df, value_name = "spearman_r")

  list(plot = plot, cor_mat = cor_mat, table = table, plot_df = cor_df, lake_data = lake_df)
}

build_abiotic_correlations <- function(
  out_dir = "correlations/outputs/abiotic",
  stats_dir = out_dir,
  data_rds = DATA_RDS_DEFAULT,
  filename_base = "abiotic_correlations"
) {
  result <- prepare_abiotic_correlation_result(data_rds = data_rds)
  files <- save_correlation_outputs(
    plot = result$plot,
    table = result$table,
    out_dir = out_dir,
    stats_dir = stats_dir,
    filename_base = filename_base,
    table_filename = "abiotic_correlation_table.csv",
    width = 9,
    height = 8
  )
  result$files <- files
  invisible(result)
}

build_all_correlations <- function(
  data_rds = DATA_RDS_DEFAULT,
  outputs_dir = "correlations/outputs",
  readme_path = file.path(outputs_dir, "README.md")
) {
  biotic <- build_biotic_correlations(
    out_dir = file.path(outputs_dir, "biotic"),
    stats_dir = file.path(outputs_dir, "biotic"),
    data_rds = data_rds
  )

  abiotic <- build_abiotic_correlations(
    out_dir = file.path(outputs_dir, "abiotic"),
    stats_dir = file.path(outputs_dir, "abiotic"),
    data_rds = data_rds
  )

  write_correlations_results_readme(
    out_path = readme_path,
    biotic_png = "biotic/biotic_correlations.png",
    biotic_table = "biotic/biotic_correlation_table.csv",
    abiotic_png = "abiotic/abiotic_correlations.png",
    abiotic_table = "abiotic/abiotic_correlation_table.csv"
  )

  invisible(list(
    biotic = biotic,
    abiotic = abiotic,
    readme = readme_path
  ))
}
