## Model 16 subset: deforestation -> eutrophication -> cyclopoids
##
## Purpose:
## - reduced version of model 16 keeping only deforestation, the eutrophication
##   composite, and cyclopoids
## - includes both a direct deforestation -> cyclopoids path and an indirect
##   deforestation -> eutrophication -> cyclopoids path

MODEL_NAME <- "model16_subset_direct_eutroph"

try({
  if (exists("df_m")) {
    z <- function(x) as.numeric(scale(as.numeric(x)))
    comp <- cbind(
      chl = z(df_m$chl_log1p),
      temp = z(df_m$temp),
      oxy_rev = -z(df_m$oxy),
      pH_rev = -z(df_m$pH)
    )
    df_m$eutroph <- as.numeric(scale(rowMeans(comp)))
  }
}, silent = TRUE)

bf_eut <- brms::bf(eutroph ~ defor, family = stats::gaussian())
bf_cyc <- brms::bf(cyc_log1p ~ defor + eutroph, family = stats::gaussian())

PATHS <- list(
  c("defor", "eutroph"),
  c("eutroph", "cyclog1p"),
  c("defor", "cyclog1p")
)

EFFECT_PANELS <- c(
  "defor->eutroph",
  "eutroph->cyc_log1p",
  "defor->cyc_log1p"
)

pri <- c(
  brms::prior(normal(0, 1), class = "b", resp = "eutroph"),
  brms::prior(normal(0, 1), class = "b", resp = "cyclog1p"),
  brms::prior(normal(0, 2), class = "Intercept", resp = "eutroph"),
  brms::prior(normal(0, 2), class = "Intercept", resp = "cyclog1p")
)

rescor_flag <- FALSE
