# SEM model: model_run



Output directory: `sem_runner/results_piecewise/model15_all_lakes_temp_do_chl_ph_cal_cyc_norescor`



## Lake coverage



Model uses 39 lake-level medians (deep/single habitat).

Deep/single candidates prior to dropping missing covariates: 39.

Unique lakes in the main dataset before exclusions: 41.



### Piecewise per-equation sample sizes



| equation | response | n_lakes |
|---|---|---|
| bf_cal | cal_log1p | 39 |
| bf_cyc | cyc_log1p | 39 |
| bf_oxy | oxy | 39 |
| bf_pH | pH | 39 |
| bf_chl0 | chl_log1p | 28 |
| bf_temp | temp | 28 |



### Missing data by submodel (piecewise)



Each submodel is fit on the subset of lakes with complete data for the response and its predictors. The table below summarizes sample sizes and the most common missing variables; per-submodel dropped-lake lists are written as `missingness_<response>.csv` in the output directory.



| equation | response | n_total | n_used | n_dropped | top_missing |
|---|---|---|---|---|---|
| bf_chl0 | chl_log1p | 39 | 28 | 11 | deforestation (% forest loss) (11); chlorophyll (0) |
| bf_temp | temp | 39 | 28 | 11 | deforestation (% forest loss) (11); temperature (0) |
| bf_cal | cal_log1p | 39 | 39 |  0 | calanoids (0); oxygen (0); pH (0) |
| bf_cyc | cyc_log1p | 39 | 39 |  0 | cyclopoids (0); pH (0); calanoids (0) |
| bf_oxy | oxy | 39 | 39 |  0 | oxygen (0); temperature (0); chlorophyll (0) |
| bf_pH | pH | 39 | 39 |  0 | pH (0); chlorophyll (0) |



## Model validation (MCMC + posterior predictive checks)



For each piecewise submodel we record basic MCMC convergence diagnostics (R-hat and effective sample sizes) and run posterior predictive checks (density overlays) plus posterior predictive PIT checks (histogram + ECDF vs Uniform(0,1)). PIT here is computed from posterior predictive draws on the fitted data (not leave-one-out).



| response | n_obs | max_rhat | min_ess_bulk | min_ess_tail | n_divergent | n_treedepth_saturated |
|---|---|---|---|---|---|---|
| cal_log1p | 39 | 1.000445 | 5130.758 | 6913.109 | 0 | 0 |
| chl_log1p | 28 | 1.001248 | 4832.511 | 6200.458 | 0 | 0 |
| cyc_log1p | 39 | 1.001362 | 4317.833 | 5600.595 | 0 | 0 |
| oxy | 39 | 1.001008 | 4969.090 | 7045.634 | 0 | 0 |
| pH | 39 | 1.000391 | 4586.129 | 6078.045 | 0 | 0 |
| temp | 28 | 1.000776 | 5050.046 | 6733.414 | 0 | 0 |



Artifacts (per response) are written to the output directory with prefixes `ppc_` (density PPC) and `ppc_pit_` (PIT histogram/ECDF), and a full parameter-level diagnostics table `mcmc_diagnostics_params_<response>.csv`.



## Composite figure



![](composite_sem_AA_BC.png)



Figure. Piecewise structural equation model (SEM) of deep/single lakes. The path diagram in panel A shows directed links among responses; edge labels are posterior median regression coefficients from BRMS Gaussian components (density variables on the natural-log scale, ln(1 + x)). Residual correlations were set to zero (rescor = FALSE).



Panel A (Path diagram). Directed graph with edges labeled by posterior median standardized coefficients (global SD scaling). Raw and standardized edge summaries are saved in `path_diagram_edges.csv` (columns `*_raw` and `*_std`). Node names: thermocline depth = “thermo”, chlorophyll = “chl_log1p”, calanoids = “cal_log1p”, cyclopoids = “cyc_log1p”, pH = “pH”, oxygen = “oxy”.



Panels B–C (conditional effects). Fitted mean responses with 95% credible bands for exemplar paths; x‑axes cover the central 96% of observed predictor values (2nd–98th percentiles).



Data and preprocessing. All abiotic predictors are water‑column integrated summaries (integrated_chl, integrated_fDOM, integrated_pH, integrated_temp, integrated_DO_percent; thermocline_depth_m from profiles). Density responses are modeled on the natural-log scale using ln(1 + x) transforms. 



## Coefficients



### coefficients_callog1p_model15_all_lakes_temp_do_chl_ph_cal_cyc_norescor.csv

| Parameter | Term | Mean | SD | Q2.5 | Q97.5 |
|---|---|---|---|---|---|
| b_Intercept | Intercept |  0.353091397 | 1.222561697 | -2.02440488 | 2.75340369 |
| b_oxy | oxygen |  0.005240471 | 0.009320022 | -0.01312064 | 0.02350267 |
| b_pH | pH | -0.037036861 | 0.191239092 | -0.41792454 | 0.33792875 |





### coefficients_chllog1p_model15_all_lakes_temp_do_chl_ph_cal_cyc_norescor.csv

| Parameter | Term | Mean | SD | Q2.5 | Q97.5 |
|---|---|---|---|---|---|
| b_Intercept | Intercept | 0.308476104 | 0.074191276 | 0.1608920182 | 0.454179096 |
| b_defor | deforestation (% forest loss) | 0.004815512 | 0.002338301 | 0.0001293036 | 0.009467657 |





### coefficients_cyclog1p_model15_all_lakes_temp_do_chl_ph_cal_cyc_norescor.csv

| Parameter | Term | Mean | SD | Q2.5 | Q97.5 |
|---|---|---|---|---|---|
| b_Intercept | Intercept |  1.34329406 | 1.4503941 | -1.4704885 | 4.1693045 |
| b_pH | pH | -0.06152863 | 0.2049609 | -0.4585138 | 0.3421109 |
| b_cal_log1p | calanoids |  0.03054826 | 0.2088654 | -0.3821353 | 0.4326469 |





### coefficients_oxy_model15_all_lakes_temp_do_chl_ph_cal_cyc_norescor.csv

| Parameter | Term | Mean | SD | Q2.5 | Q97.5 |
|---|---|---|---|---|---|
| b_Intercept | Intercept |  97.267771 | 9.5464718 |  78.3578852 | 116.085298 |
| b_temp | temperature |   1.591283 | 0.6781131 |   0.2746683 |   2.929122 |
| b_chl_log1p | chlorophyll | -24.800007 | 7.6435884 | -39.8110556 |  -9.640362 |





### coefficients_pH_model15_all_lakes_temp_do_chl_ph_cal_cyc_norescor.csv

| Parameter | Term | Mean | SD | Q2.5 | Q97.5 |
|---|---|---|---|---|---|
| b_Intercept | Intercept |  7.622735 | 0.2044147 |  7.220043 |  8.0228791 |
| b_chl_log1p | chlorophyll | -1.240814 | 0.4043533 | -2.032145 | -0.4440639 |





### coefficients_temp_model15_all_lakes_temp_do_chl_ph_cal_cyc_norescor.csv

| Parameter | Term | Mean | SD | Q2.5 | Q97.5 |
|---|---|---|---|---|---|
| b_Intercept | Intercept | 13.08589616 | 1.05249976 | 10.9969452 | 15.16793490 |
| b_defor | deforestation (% forest loss) | -0.05075343 | 0.03314741 | -0.1154758 |  0.01461298 |





## Effects (direct, indirect, total)



Effects are computed from posterior draws. Direct effects are regression slopes for each arrow. Indirect effects are computed draw-by-draw as products of coefficients along each directed path; total effects are direct + summed indirect (when applicable).



### Total effects (`path_effects_total.csv`)

| direct_included | n_indirect_paths | median | q2.5 | q97.5 | n_draws | nonzero_95 | median_std | q2.5_std | q97.5_std | nonzero_95_std | from_label | to_label |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| TRUE | 0 |  3.308764e-02 | -3.821353e-01 |  4.326469e-01 | 12000 | FALSE |  0.0276862261 | -0.319753403 |  0.362019218 | FALSE | calanoids | cyclopoids |
| FALSE | 2 | -8.293317e-02 | -6.198414e-01 |  4.372300e-01 | 12000 | FALSE | -0.0312898743 | -0.233860098 |  0.164962583 | FALSE | chlorophyll | calanoids |
| FALSE | 3 |  6.469888e-02 | -4.632031e-01 |  6.425019e-01 | 12000 | FALSE |  0.0204253880 | -0.146232872 |  0.202837355 | FALSE | chlorophyll | cyclopoids |
| TRUE | 0 | -2.472773e+01 | -3.981106e+01 | -9.640362e+00 | 12000 | TRUE | -0.4423122750 | -0.712112315 | -0.172440061 | TRUE | chlorophyll | oxygen |
| TRUE | 0 | -1.238944e+00 | -2.032145e+00 | -4.440639e-01 | 12000 | TRUE | -0.4586235752 | -0.752245031 | -0.164380440 | TRUE | chlorophyll | pH |
| FALSE | 3 | -6.563091e-04 | -5.202463e-03 |  2.875347e-03 | 12000 | FALSE | -0.0167115957 | -0.132470289 |  0.073214933 | FALSE | deforestation | calanoids |
| TRUE | 0 |  4.811806e-03 |  1.293036e-04 |  9.467657e-03 | 12000 | TRUE |  0.3247446744 |  0.008726588 |  0.638964081 | TRUE | deforestation | chlorophyll |
| FALSE | 4 |  2.107559e-04 | -2.573315e-03 |  3.762673e-03 | 12000 | FALSE |  0.0044904230 | -0.054827755 |  0.080168537 | FALSE | deforestation | cyclopoids |
| FALSE | 2 | -1.929343e-01 | -4.069233e-01 | -3.138705e-02 | 12000 | TRUE | -0.2329100839 | -0.491237305 | -0.037890404 | TRUE | deforestation | oxygen |
| FALSE | 1 | -5.555584e-03 | -1.421622e-02 | -6.927397e-05 | 12000 | TRUE | -0.1387932478 | -0.355158969 | -0.001730648 | TRUE | deforestation | pH |
| TRUE | 0 | -5.082655e-02 | -1.154758e-01 |  1.461298e-02 | 12000 | FALSE | -0.2994865068 | -0.680420806 |  0.086104436 | FALSE | deforestation | temperature |
| TRUE | 0 |  5.328995e-03 | -1.312064e-02 |  2.350267e-02 | 12000 | FALSE |  0.1124025257 | -0.276748822 |  0.495733174 | FALSE | oxygen | calanoids |
| FALSE | 1 |  3.992273e-05 | -4.593268e-03 |  5.322496e-03 | 12000 | FALSE |  0.0007046102 | -0.081068191 |  0.093938593 | FALSE | oxygen | cyclopoids |
| TRUE | 0 | -3.369538e-02 | -4.179245e-01 |  3.379287e-01 | 12000 | FALSE | -0.0343432331 | -0.425959922 |  0.344426063 | FALSE | pH | calanoids |
| TRUE | 1 | -6.336247e-02 | -4.706480e-01 |  3.508801e-01 | 12000 | FALSE | -0.0540382097 | -0.401388671 |  0.299245464 | FALSE | pH | cyclopoids |
| FALSE | 1 |  6.806668e-03 | -2.236656e-02 |  4.453580e-02 | 12000 | FALSE |  0.0294142224 | -0.096654477 |  0.192456269 | FALSE | temperature | calanoids |
| FALSE | 1 |  4.827971e-05 | -7.728588e-03 |  9.097671e-03 | 12000 | FALSE |  0.0001745764 | -0.027946085 |  0.032896600 | FALSE | temperature | cyclopoids |
| TRUE | 0 |  1.589040e+00 |  2.746683e-01 |  2.929122e+00 | 12000 | TRUE |  0.3255568447 |  0.056273069 |  0.600108149 | TRUE | temperature | oxygen |



### Direct effects (`path_effects.csv`)

| path | median | q2.5 | q97.5 | n_draws | nonzero_95 | median_std | q2.5_std | q97.5_std | nonzero_95_std |
|---|---|---|---|---|---|---|---|---|---|
| deforestation -> temperature |  -0.050826546 | -1.154758e-01 |  0.014612982 | 12000 | FALSE | -0.29948651 | -0.680420806 |  0.08610444 | FALSE |
| deforestation -> chlorophyll |   0.004811806 |  1.293036e-04 |  0.009467657 | 12000 | TRUE |  0.32474467 |  0.008726588 |  0.63896408 | TRUE |
| temperature -> oxygen |   1.589039740 |  2.746683e-01 |  2.929121944 | 12000 | TRUE |  0.32555684 |  0.056273069 |  0.60010815 | TRUE |
| chlorophyll -> oxygen | -24.727726508 | -3.981106e+01 | -9.640362493 | 12000 | TRUE | -0.44231227 | -0.712112315 | -0.17244006 | TRUE |
| chlorophyll -> pH |  -1.238944107 | -2.032145e+00 | -0.444063909 | 12000 | TRUE | -0.45862358 | -0.752245031 | -0.16438044 | TRUE |
| oxygen -> calanoids |   0.005328995 | -1.312064e-02 |  0.023502673 | 12000 | FALSE |  0.11240253 | -0.276748822 |  0.49573317 | FALSE |
| pH -> calanoids |  -0.033695376 | -4.179245e-01 |  0.337928748 | 12000 | FALSE | -0.03434323 | -0.425959922 |  0.34442606 | FALSE |
| pH -> cyclopoids |  -0.062531360 | -4.585138e-01 |  0.342110931 | 12000 | FALSE | -0.05332941 | -0.391040132 |  0.29176677 | FALSE |
| calanoids -> cyclopoids |   0.033087637 | -3.821353e-01 |  0.432646926 | 12000 | FALSE |  0.02768623 | -0.319753403 |  0.36201922 | FALSE |



### Indirect effects (`path_effects_indirect.csv`)



Indirect effects are reported per directed path (with a `via` and explicit `path`). See `path_effects_indirect.csv` for the full table.





## Figures



### Posterior predictive checks

![](ppc_cal_log1p.png)

![](ppc_chl_log1p.png)

![](ppc_cyc_log1p.png)

![](ppc_oxy.png)

![](ppc_pH.png)

![](ppc_pit_ecdf_cal_log1p.png)

![](ppc_pit_ecdf_chl_log1p.png)

![](ppc_pit_ecdf_cyc_log1p.png)

![](ppc_pit_ecdf_oxy.png)

![](ppc_pit_ecdf_pH.png)

![](ppc_pit_ecdf_temp.png)

![](ppc_pit_hist_cal_log1p.png)

![](ppc_pit_hist_chl_log1p.png)

![](ppc_pit_hist_cyc_log1p.png)

![](ppc_pit_hist_oxy.png)

![](ppc_pit_hist_pH.png)

![](ppc_pit_hist_temp.png)

![](ppc_temp.png)



### Path diagram

![](path_diagram.png)



