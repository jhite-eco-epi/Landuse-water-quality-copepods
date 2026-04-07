## Model 16 subset: deforestation -> cyclopoids only
##
## Purpose:
## - minimal comparison model with a single direct path

MODEL_NAME <- "model16_subset_direct_only"

bf_cyc <- brms::bf(cyc_log1p ~ defor, family = stats::gaussian())

PATHS <- list(
  c("defor", "cyclog1p")
)

EFFECT_PANELS <- c(
  "defor->cyc_log1p"
)

pri <- c(
  brms::prior(normal(0, 1), class = "b", resp = "cyclog1p"),
  brms::prior(normal(0, 2), class = "Intercept", resp = "cyclog1p")
)

rescor_flag <- FALSE
