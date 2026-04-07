# SEM model: model_run



Output directory: `sem_runner/results_piecewise/model16_subset_direct_eutroph`



## Lake coverage



Model uses 28 lake-level medians (deep/single habitat).

Deep/single candidates prior to dropping missing covariates: 39.

Unique lakes in the main dataset before exclusions: 41.



### Piecewise per-equation sample sizes



| equation | response | n_lakes |
|---|---|---|
| bf_cyc | cyc_log1p | 28 |
| bf_eut | eutroph | 28 |



### Missing data by submodel (piecewise)



Each submodel is fit on the subset of lakes with complete data for the response and its predictors. The table below summarizes sample sizes and the most common missing variables; per-submodel dropped-lake lists are written as `missingness_<response>.csv` in the output directory.



| equation | response | n_total | n_used | n_dropped | top_missing |
|---|---|---|---|---|---|
| bf_cyc | cyc_log1p | 39 | 28 | 11 | deforestation (% forest loss) (11); cyclopoids (0); eutrophication composite (0) |
| bf_eut | eutroph | 39 | 28 | 11 | deforestation (% forest loss) (11); eutrophication composite (0) |



## Model validation (MCMC + posterior predictive checks)



For each piecewise submodel we record basic MCMC convergence diagnostics (R-hat and effective sample sizes) and run posterior predictive checks (density overlays) plus posterior predictive PIT checks (histogram + ECDF vs Uniform(0,1)). PIT here is computed from posterior predictive draws on the fitted data (not leave-one-out).



| response | n_obs | max_rhat | min_ess_bulk | min_ess_tail | n_divergent | n_treedepth_saturated |
|---|---|---|---|---|---|---|
| cyc_log1p | 28 | 1.001563 | 4169.353 | 5776.905 | 0 | 0 |
| eutroph | 28 | 1.000667 | 4621.502 | 5327.602 | 0 | 0 |



Artifacts (per response) are written to the output directory with prefixes `ppc_` (density PPC) and `ppc_pit_` (PIT histogram/ECDF), and a full parameter-level diagnostics table `mcmc_diagnostics_params_<response>.csv`.



## Composite figure



![](composite_sem_AA_BC.png)



Figure. Piecewise structural equation model (SEM) of deep/single lakes. The path diagram in panel A shows directed links among responses; edge labels are posterior median regression coefficients from BRMS Gaussian components (density variables on the natural-log scale, ln(1 + x)). Residual correlations were set to zero (rescor = FALSE).



Panel A (Path diagram). Directed graph with edges labeled by posterior median standardized coefficients (global SD scaling). Raw and standardized edge summaries are saved in `path_diagram_edges.csv` (columns `*_raw` and `*_std`). Node names: thermocline depth = “thermo”, chlorophyll = “chl_log1p”, calanoids = “cal_log1p”, cyclopoids = “cyc_log1p”, pH = “pH”, oxygen = “oxy”.



Panels B–C (conditional effects). Fitted mean responses with 95% credible bands for exemplar paths; x‑axes cover the central 96% of observed predictor values (2nd–98th percentiles).



Data and preprocessing. All abiotic predictors are water‑column integrated summaries (integrated_chl, integrated_fDOM, integrated_pH, integrated_temp, integrated_DO_percent; thermocline_depth_m from profiles). Density responses are modeled on the natural-log scale using ln(1 + x) transforms. 



## Coefficients



### coefficients_cyclog1p_model16_subset_direct_eutroph.csv

| Parameter | Term | Mean | SD | Q2.5 | Q97.5 |
|---|---|---|---|---|---|
| b_Intercept | Intercept |  0.49912332 | 0.35989351 | -0.200519790 | 1.21373046 |
| b_defor | deforestation (% forest loss) |  0.01790746 | 0.01163709 | -0.005259017 | 0.04047795 |
| b_eutroph | eutrophication composite | -0.03346644 | 0.22074852 | -0.469075829 | 0.40574807 |





### coefficients_eutroph_model16_subset_direct_eutroph.csv

| Parameter | Term | Mean | SD | Q2.5 | Q97.5 |
|---|---|---|---|---|---|
| b_Intercept | Intercept | -0.95495430 | 0.275208058 | -1.49272692 | -0.41149417 |
| b_defor | deforestation (% forest loss) |  0.03280497 | 0.008570162 |  0.01584891 |  0.04970369 |





## Effects (direct, indirect, total)



Effects are computed from posterior draws. Direct effects are regression slopes for each arrow. Indirect effects are computed draw-by-draw as products of coefficients along each directed path; total effects are direct + summed indirect (when applicable).



### Total effects (`path_effects_total.csv`)

| direct_included | n_indirect_paths | median | q2.5 | q97.5 | n_draws | nonzero_95 | median_std | q2.5_std | q97.5_std | nonzero_95_std | from_label | to_label |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| TRUE | 1 |  0.01683545 | -0.001649495 | 0.03511246 | 12000 | FALSE |  0.35870072 | -0.03514458 | 0.7481156 | FALSE | deforestation | cyclopoids |
| TRUE | 0 |  0.03276324 |  0.015848913 | 0.04970369 | 12000 | TRUE |  0.45044741 |  0.21789973 | 0.6833542 | TRUE | deforestation | eutrophication composite |
| TRUE | 0 | -0.03586464 | -0.469075829 | 0.40574807 | 12000 | FALSE | -0.05557978 | -0.72693122 | 0.6287916 | FALSE | eutrophication composite | cyclopoids |



### Direct effects (`path_effects.csv`)

| path | median | q2.5 | q97.5 | n_draws | nonzero_95 | median_std | q2.5_std | q97.5_std | nonzero_95_std |
|---|---|---|---|---|---|---|---|---|---|
| deforestation -> eutrophication composite |  0.03276324 |  0.015848913 | 0.04970369 | 12000 | TRUE |  0.45044741 |  0.2178997 | 0.6833542 | TRUE |
| eutrophication composite -> cyclopoids | -0.03586464 | -0.469075829 | 0.40574807 | 12000 | FALSE | -0.05557978 | -0.7269312 | 0.6287916 | FALSE |
| deforestation -> cyclopoids |  0.01800832 | -0.005259017 | 0.04047795 | 12000 | FALSE |  0.38369011 | -0.1120500 | 0.8624341 | FALSE |



### Indirect effects (`path_effects_indirect.csv`)



Indirect effects are reported per directed path (with a `via` and explicit `path`). See `path_effects_indirect.csv` for the full table.





## Figures



### Posterior predictive checks

![](ppc_cyc_log1p.png)

![](ppc_eutroph.png)

![](ppc_pit_ecdf_cyc_log1p.png)

![](ppc_pit_ecdf_eutroph.png)

![](ppc_pit_hist_cyc_log1p.png)

![](ppc_pit_hist_eutroph.png)



### Path diagram

![](path_diagram.png)



