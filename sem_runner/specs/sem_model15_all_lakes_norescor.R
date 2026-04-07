## Model 15 (all lakes): same as model 15, but intended to be run
## with --keep-all-lakes=TRUE so no a priori lake exclusions are applied.

MODEL_NAME <- "model15_all_lakes_temp_do_chl_ph_cal_cyc_norescor"

# deforestation -> temperature
bf_temp <- brms::bf(temp ~ defor, family = stats::gaussian())

# deforestation -> chlorophyll
bf_chl0 <- brms::bf(chl_log1p ~ defor, family = stats::gaussian())

# temperature + chlorophyll -> oxygen
bf_oxy  <- brms::bf(oxy ~ temp + chl_log1p, family = stats::gaussian())

# chlorophyll -> pH
bf_pH   <- brms::bf(pH ~ chl_log1p, family = stats::gaussian())

# oxygen + pH -> calanoids
bf_cal  <- brms::bf(cal_log1p ~ oxy + pH, family = stats::gaussian())

# pH + calanoids -> cyclopoids
bf_cyc  <- brms::bf(cyc_log1p ~ pH + cal_log1p, family = stats::gaussian())

bf_chl <- bf_temp + bf_chl0 + bf_oxy + bf_pH + bf_cal
bf_cyc <- bf_cyc

PATHS <- list(
  c("defor","temp"),
  c("defor","chllog1p"),
  c("temp","oxy"),
  c("chllog1p","oxy"),
  c("chllog1p","pH"),
  c("oxy","callog1p"),
  c("pH","callog1p"),
  c("pH","cyclog1p"),
  c("callog1p","cyclog1p")
)

EFFECT_PANELS <- c(
  "defor->temp",
  "defor->chl_log1p",
  "temp->oxy",
  "chl_log1p->oxy",
  "chl_log1p->pH",
  "oxy->cal_log1p",
  "pH->cal_log1p",
  "pH->cyc_log1p",
  "cal_log1p->cyc_log1p"
)

pri <- c(
  brms::prior(normal(0, 1), class = "b", resp = "chllog1p"),
  brms::prior(normal(0, 1), class = "b", resp = "callog1p"),
  brms::prior(normal(0, 1), class = "b", resp = "cyclog1p"),
  brms::prior(normal(0, 1), class = "b", resp = "pH"),
  brms::prior(normal(0, 1), class = "b", resp = "oxy"),
  brms::prior(normal(0, 1), class = "b", resp = "temp")
)

rescor_flag <- FALSE
