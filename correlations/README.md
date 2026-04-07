# Correlations

Reusable scripts for generating lake-level correlation heatmaps for:

- `biotic`: zooplankton taxa densities across lakes
- `abiotic`: integrated environmental predictors across lakes

Both workflows share code in `correlations/helpers_correlations.R`.

## Run

From the project root:

```bash
Rscript correlations/run_all_correlations.R
```

The combined runner regenerates both figures and writes a results-focused markdown report to `correlations/outputs/README.md`.

To execute individually:

```bash
Rscript correlations/make_biotic_correlations.R
```
and / or:

 ```bash
Rscript correlations/make_abiotic_correlations.R
 ```

## Outputs

- `correlations/outputs/biotic/biotic_correlations.png`
- `correlations/outputs/biotic/biotic_correlations.tiff`
- `correlations/outputs/biotic/biotic_correlation_table.csv`
- `correlations/outputs/abiotic/abiotic_correlations.png`
- `correlations/outputs/abiotic/abiotic_correlations.tiff`
- `correlations/outputs/abiotic/abiotic_correlation_table.csv`
- `correlations/outputs/README.md`
