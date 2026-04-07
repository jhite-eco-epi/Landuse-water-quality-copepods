# Data Pipeline (Sample-Level Main)

This folder contains the active, sample-level data integration pipeline used to build the main dataset for all downstream analyses.

## What it does
- Reads raw inputs (sonde profiles, zooplankton counts, bathymetry, geography)
- Detects sonde deployments and classifies shallow/deep per lake
- Computes depth-weighted (trapezoidal) deployment means for chemistry (integrated_*)
- Converts zooplankton counts to densities (individuals per liter) per sample
- Joins lake-level variables and spatial attributes
- Matches zooplankton samples to sonde deployments (with lake-average fallback)
- Writes the sample-level main dataset

## Inputs
- `data/Sonde Data/VancouverSondeData2023_UTF8.csv`
- `data/Bathymetry Data/LakeBathymetry.csv` (regenerated from Excel)
- `data/LakeGeography/Field_2023_Lake_Level_Ocean_Distance.csv`
- `data/Zooplankton ID Data/ZoopIDDataUpdated.csv`
- `data/Infection Data/42_Lake_Copepod_Infection.csv` (joined by lake + shallow/deep; single samples fall back to deep then shallow)

## How to run
```r
# From project root
Rscript run_pipeline.R
```

## Outputs (main)
- `data/main_zooplankton_data.csv`
- `data/main_zooplankton_data.rds`

## Key columns
- IDs/meta: `SampleID`, `LakeID`, `Depth_Tow_m`, `zoop_sample_type` (shallow/deep/single)
- Zooplankton: `<Taxon>_density_per_L` (per-sample densities)
- Sonde (depth-weighted): `integrated_temp`, `integrated_pH`, `integrated_DO_percent`, `integrated_DO_mgL`, `integrated_chl`, `integrated_phyco`, `integrated_fDOM`, plus `surface_*`, `bottom_*`
- Lake-level: `MaxDepth`, `MeanDepth`, `SurfaceArea_m2`, `Volume_m3` (if available), `depth_ratio`
- Geography: `Latitude`, `Longitude`, `DFO_m`, `Lake_Level_m`
- Infection: `infection_*` (e.g., `infection_burden`, `infection_cope_counts`, `infection_worm_counts`)
- Matching quality: `sonde_match_quality` (deployment_matched / lake_average / no_sonde_data)

## Notes
- No lake averaging is performed here; sample-level is the authoritative main. Any averaging or selection (e.g., deep-only) should be done in the analysis stage.
- The legacy lake-averaged pipeline has been deprecated and removed.
