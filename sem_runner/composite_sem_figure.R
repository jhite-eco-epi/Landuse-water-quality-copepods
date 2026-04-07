#!/usr/bin/env Rscript
# Composite SEM figure using ggplot objects sourced from helper scripts
# Layout:
#   A A
#   B C

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# Allow overriding which effect panels to show
SEM_EFFECT_B <- Sys.getenv("SEM_EFFECT_B", unset = "chl_log1p->cyc_log1p") # predictor->response
SEM_EFFECT_C <- Sys.getenv("SEM_EFFECT_C", unset = "pH->chl_log1p")
OUT_DIR <- Sys.getenv("SEM_OUT_DIR", unset = "sem_runner/results_deep")

# If a spec is provided, let it define EFFECT_PANELS to drive effect panels
SPEC_ENV <- Sys.getenv("SEM_SPEC_PATH", unset = "")
MODEL_NAME <- Sys.getenv("SEM_MODEL_NAME", unset = "")
if (nzchar(SPEC_ENV) && file.exists(SPEC_ENV)) {
  try(sys.source(SPEC_ENV, envir = environment()), silent = TRUE)
} else if (nzchar(MODEL_NAME)) {
  cand <- file.path("sem_runner", "specs", paste0("sem_", MODEL_NAME, ".R"))
  if (file.exists(cand)) try(sys.source(cand, envir = environment()), silent = TRUE)
}
if (exists("EFFECT_PANELS", inherits = FALSE) && is.character(EFFECT_PANELS)) {
  EFFECT_LIST <- EFFECT_PANELS
} else {
  EFFECT_LIST <- c(SEM_EFFECT_B, SEM_EFFECT_C)
}

# If we still have fewer than 3 panels and PATHS exist, derive effects from PATHS
if (length(EFFECT_LIST) < 3 && exists("PATHS", inherits = FALSE) && length(PATHS) > 0) {
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
      from_node <- pth[k]
      to_node   <- pth[k + 1]
      derived <- c(derived, paste0(term_of(from_node), "->", term_of(to_node)))
    }
  }
  derived <- unique(derived)
  if (length(derived) > length(EFFECT_LIST)) {
    EFFECT_LIST <- derived
  }
}

message("SEM composite effects (B..): ", paste(EFFECT_LIST, collapse = ", "))

# Source helper plot objects
source("sem_runner/scripts/plot_sem_path_panel.R")   # exports p_path
source("sem_runner/scripts/plot_sem_effect_panel.R") # exports function make_effect_plot(effect_spec)

# Build panels
pA <- p_path
p_effects <- lapply(EFFECT_LIST, make_effect_plot)

# Add panel tags aligned to inner axes (top-left)
panel_label <- function(p, lbl) {
  p + annotate("text", x = -Inf, y = Inf, label = lbl,
               hjust = -0.15, vjust = 1.25, size = 7, fontface = "bold")
}
pA <- panel_label(pA, "A")
if (length(p_effects) > 0) {
  letters <- LETTERS[2:(1 + length(p_effects))]
  for (i in seq_along(p_effects)) {
    p_effects[[i]] <- panel_label(p_effects[[i]], letters[i])
  }
}

# Layout: A full width at top; effects stacked in 2 columns below
rows <- ceiling(length(p_effects) / 2)
if (length(p_effects) == 0) {
  final_plot <- pA & theme(plot.background = element_blank())
} else {
  effects_grid <- wrap_plots(p_effects, ncol = 2)
  final_plot <- (pA + theme(plot.margin = margin(0,0,6,0))) /
    (effects_grid + theme(plot.margin = margin(0,0,0,0)))
  final_plot <- final_plot & theme(plot.background = element_blank(), plot.margin = margin(2,2,2,2))
}

print(final_plot)

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

base_h <- 6
per_row <- 3
plot_h <- if (length(p_effects) == 0) base_h else base_h + rows * per_row
ggsave(file.path(OUT_DIR, "composite_sem_AA_BC.tiff"),
       final_plot, width = 14, height = plot_h, dpi = 600, bg = "white", compression = "lzw")
ggsave(file.path(OUT_DIR, "composite_sem_AA_BC.png"),
       final_plot, width = 14, height = plot_h, dpi = 175, bg = "white")

message("Wrote: ", file.path(OUT_DIR, "composite_sem_AA_BC.(tiff|png)"))


