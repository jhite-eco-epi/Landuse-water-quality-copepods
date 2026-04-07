# Data Pipeline for Zooplankton Density Analysis

# Load required libraries
library(tidyverse)
library(here)

# Source the sample-level data processing functions
source("data_pipeline/data_processing_functions.R")

# Set data paths
data_paths <- list(
  sonde = "01_WorkingData/Sonde Data/VancouverSondeData2023_UTF8.csv",
  bathy = "01_WorkingData/Bathymetry Data/LakeBathymetry.csv",
  geo = "01_WorkingData/LakeGeography/Field_2023_Lake_Level_Ocean_Distance.csv",
  zoop = "01_WorkingData/Zooplankton ID Data/ZoopIDDataUpdated.csv",
  infection = "01_WorkingData/Infection Data/42_Lake_Copepod_Infection.csv"
)

# Run the pipeline
message("Starting sample-level data pipeline...")
message("====================================")

sample_df <- create_sample_level_dataframe(
  sonde_path = data_paths$sonde,
  bathy_path = data_paths$bathy,
  geo_path = data_paths$geo,
  zoop_path = data_paths$zoop,
  infection_path = data_paths$infection
)

# Save as main dataset (sample-level is our main)
output_path <- "01_WorkingData/main_zooplankton_data.csv"
write.csv(sample_df, output_path, row.names = FALSE)
message("\nMain dataframe saved to: ", output_path)

# Also save as RDS for faster loading in R
rds_path <- "01_WorkingData/main_zooplankton_data.rds"
saveRDS(sample_df, rds_path)
message("RDS file saved to: ", rds_path)

message("\nPipeline complete!")
