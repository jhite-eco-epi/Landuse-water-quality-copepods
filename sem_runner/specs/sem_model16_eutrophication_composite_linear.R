## Model 16: Deforestation -> eutrophication composite -> calanoids / cyclopoids
##
## Purpose:
## - collapse correlated water-quality mediators into one standardized
##   eutrophication composite
## - retain calanoids as a biotic mediator
## - add a direct deforestation -> cyclopoids path
##
## Composite orientation:
## - higher values indicate more eutrophic conditions
## - positive loadings: chlorophyll, temperature
## - negative loadings: oxygen, pH

MODEL_NAME <- "model16_eutrophication_composite_linear"

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
bf_cal <- brms::bf(cal_log1p ~ eutroph, family = stats::gaussian())
bf_cyc <- brms::bf(cyc_log1p ~ eutroph + cal_log1p + defor, family = stats::gaussian())

bf_chl <- bf_eut + bf_cal
bf_cyc <- bf_cyc

PATHS <- list(
  c("defor", "eutroph"),
  c("eutroph", "callog1p"),
  c("eutroph", "cyclog1p"),
  c("callog1p", "cyclog1p"),
  c("defor", "cyclog1p")
)

EFFECT_PANELS <- c(
  "defor->eutroph",
  "eutroph->cal_log1p",
  "eutroph->cyc_log1p",
  "cal_log1p->cyc_log1p",
  "defor->cyc_log1p"
)

pri <- c(
  brms::prior(normal(0, 1), class = "b", resp = "eutroph"),
  brms::prior(normal(0, 1), class = "b", resp = "callog1p"),
  brms::prior(normal(0, 1), class = "b", resp = "cyclog1p"),
  brms::prior(normal(0, 2), class = "Intercept", resp = "eutroph"),
  brms::prior(normal(0, 2), class = "Intercept", resp = "callog1p"),
  brms::prior(normal(0, 2), class = "Intercept", resp = "cyclog1p")
)

rescor_flag <- FALSE
