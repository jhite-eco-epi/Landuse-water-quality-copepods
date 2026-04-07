#!/usr/bin/env Rscript
# Exports a ggplot object `p_path` representing the SEM path diagram.
#
# In this project, we always run the SEM on deep/single lake medians upstream.
# The diagram PNG is rendered by `sem_runner/path_effects_and_diagram.R` and
# stored as `path_diagram.png` in the run output directory; this script simply
# embeds that PNG as a ggplot for composite figure layout.

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

OUT_DIR <- Sys.getenv("SEM_OUT_DIR", unset = "sem_runner/results_deep")
DIAGRAM_PATH <- file.path(OUT_DIR, "path_diagram.png")

if (!requireNamespace("png", quietly = TRUE)) {
  stop("Package 'png' is required to read the path diagram PNG.", call. = FALSE)
}
if (!file.exists(DIAGRAM_PATH)) {
  stop("Expected path diagram PNG not found at ", DIAGRAM_PATH, call. = FALSE)
}

img <- png::readPNG(DIAGRAM_PATH)
p_path <- ggplot() +
  annotation_custom(grid::rasterGrob(img, interpolate = TRUE), xmin = 0, xmax = 2, ymin = 0, ymax = 1) +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_fixed(ratio = 0.5, clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))
