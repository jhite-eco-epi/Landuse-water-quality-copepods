# SEM model: model_run



Output directory: `sem_runner/results_piecewise/model16_subset_direct_only`



## Lake coverage



Model uses 28 lake-level medians (deep/single habitat).

Deep/single candidates prior to dropping missing covariates: 39.

Unique lakes in the main dataset before exclusions: 41.



### Piecewise per-equation sample sizes



| equation | response | n_lakes |
|---|---|---|
| bf_cyc | cyc_log1p | 28 |



### Missing data by submodel (piecewise)



Each submodel is fit on the subset of lakes with complete data for the response and its predictors. The table below summarizes sample sizes and the most common missing variables; per-submodel dropped-lake lists are written as `missingness_<response>.csv` in the output directory.



| equation | response | n_total | n_used | n_dropped | top_missing |
|---|---|---|---|---|---|
| bf_cyc | cyc_log1p | 39 | 28 | 11 | deforestation (% forest loss) (11); cyclopoids (0) |



## Model validation (MCMC + posterior predictive checks)



For each piecewise submodel we record basic MCMC convergence diagnostics (R-hat and effective sample sizes) and run posterior predictive checks (density overlays) plus posterior predictive PIT checks (histogram + ECDF vs Uniform(0,1)). PIT here is computed from posterior predictive draws on the fitted data (not leave-one-out).



| response | n_obs | max_rhat | min_ess_bulk | min_ess_tail | n_divergent | n_treedepth_saturated |
|---|---|---|---|---|---|---|
| cyc_log1p | 28 | 1.000339 | 4502.24 | 5223.884 | 0 | 0 |



Artifacts (per response) are written to the output directory with prefixes `ppc_` (density PPC) and `ppc_pit_` (PIT histogram/ECDF), and a full parameter-level diagnostics table `mcmc_diagnostics_params_<response>.csv`.



## Composite figure



![](composite_sem_AA_BC.png)



Figure. Piecewise structural equation model (SEM) of deep/single lakes. The path diagram in panel A shows directed links among responses; edge labels are posterior median regression coefficients from BRMS Gaussian components (density variables on the natural-log scale, ln(1 + x)). Residual correlations were set to zero (rescor = FALSE).



Panel A (Path diagram). Directed graph with edges labeled by posterior median standardized coefficients (global SD scaling). Raw and standardized edge summaries are saved in `path_diagram_edges.csv` (columns `*_raw` and `*_std`). Node names: thermocline depth = “thermo”, chlorophyll = “chl_log1p”, calanoids = “cal_log1p”, cyclopoids = “cyc_log1p”, pH = “pH”, oxygen = “oxy”.



Panels B–C (conditional effects). Fitted mean responses with 95% credible bands for exemplar paths; x‑axes cover the central 96% of observed predictor values (2nd–98th percentiles).



Data and preprocessing. All abiotic predictors are water‑column integrated summaries (integrated_chl, integrated_fDOM, integrated_pH, integrated_temp, integrated_DO_percent; thermocline_depth_m from profiles). Density responses are modeled on the natural-log scale using ln(1 + x) transforms. 



## Coefficients



### coefficients_cyclog1p_model16_subset_direct_only.csv

| Parameter | Term | Mean | SD | Q2.5 | Q97.5 |
|---|---|---|---|---|---|
| b_Intercept | Intercept | 0.52082784 | 0.286884326 | -0.049616140 | 1.09052636 |
| b_defor | deforestation (% forest loss) | 0.01712096 | 0.008962225 | -0.000764362 | 0.03495322 |





## Effects (direct, indirect, total)



Effects are computed from posterior draws. Direct effects are regression slopes for each arrow. Indirect effects are computed draw-by-draw as products of coefficients along each directed path; total effects are direct + summed indirect (when applicable).



### Total effects (`path_effects_total.csv`)

| direct_included | n_indirect_paths | median | q2.5 | q97.5 | n_draws | nonzero_95 | median_std | q2.5_std | q97.5_std | nonzero_95_std | from_label | to_label |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| TRUE | 0 | 0.01707086 | -0.000764362 | 0.03495322 | 12000 | FALSE | 0.3637164 | -0.01628571 | 0.7447229 | FALSE | deforestation | cyclopoids |



### Direct effects (`path_effects.csv`)

| path | median | q2.5 | q97.5 | n_draws | nonzero_95 | median_std | q2.5_std | q97.5_std | nonzero_95_std |
|---|---|---|---|---|---|---|---|---|---|
| deforestation -> cyclopoids | 0.01707086 | -0.000764362 | 0.03495322 | 12000 | FALSE | 0.3637164 | -0.01628571 | 0.7447229 | FALSE |



### Indirect effects (`path_effects_indirect.csv`)



Indirect effects are reported per directed path (with a `via` and explicit `path`). See `path_effects_indirect.csv` for the full table.





## Figures



### Posterior predictive checks

![](ppc_cyc_log1p.png)

![](ppc_pit_ecdf_cyc_log1p.png)

![](ppc_pit_hist_cyc_log1p.png)



### Path diagram

![](path_diagram.png)



