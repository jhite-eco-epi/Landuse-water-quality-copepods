# Landuse-water-quality-copepods

This project investigates the impact of deforestation on Cyclopoid densities using field data.  

For ease of use the reository is self-contained.  Raw data and generated analyses are committed alongside the code.

There are three functional pieces included in this repository:
1.  Raw data and an associated data pipeline to create a main dataset
2.  Spearman correlation analysis
3.  Piecewise Bayesian Structure Equation Modeling using the brms package

#  Data 

Raw and processed data is found in the `data` directory.  To run the datapipleine that cleans, transforms, and joins raw data into a main dataframe:

```
Rscript run_pipeline.R
```

This will generate 

`/data/main_zooplankton_data.csv`
and 
`/data/main_zooplankton_data.rds`

which are used for all donwtream analyses.

# Correlations

The exploratory corelation analysis generatees Spearman correlations between measured variables.  The analysis is split into abiotic vs biotic.

To run the correlations:

```
Rscript correlations/run_all_correlations.R
```

Generated [correlations README](correlations/README.md)

# Piecewise Structured Equation Models

The SEM tests our hypothesis.  To run the SEM model:

```
Rscript sem_runner/sem_runner_piecewise.R --spec=sem_runner/specs/sem_model15.R
```

Generated [SEM README](sem_runner/results_piecewise/model15_logging_temp_do_chl_ph_cal_cyc_norescor/README.md)
