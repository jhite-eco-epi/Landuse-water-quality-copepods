## SEM spec template
## Define/override: bf_chl, bf_cyc, pri, rescor_flag
## This file is sourced by sem_runner_deep.R when passed via --spec=... and can
## override the default model components. Keep object names exactly as below.
##
## Example: Gaussian on log1p for both outcomes with residual correlation.

bf_chl <- brms::bf(chl_log1p ~ thermo + pH, family = stats::gaussian())
bf_cyc <- brms::bf(cyc_log1p ~ chl_log1p + thermo, family = stats::gaussian())

pri <- c(
  brms::prior(normal(0, 1), class = "b", resp = "chllog1p"),
  brms::prior(normal(0, 1), class = "b", resp = "cyclog1p"),
  brms::prior(normal(0, 2), class = "Intercept", resp = "chllog1p"),
  brms::prior(normal(0, 2), class = "Intercept", resp = "cyclog1p")
)

# Set to TRUE to estimate residual correlation in multivariate Gaussian models.
rescor_flag <- FALSE


