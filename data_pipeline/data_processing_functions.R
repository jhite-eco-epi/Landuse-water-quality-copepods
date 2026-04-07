# Sample-Level Data Processing Functions for Zooplankton Density Analysis
# These functions create a dataset with individual samples (no lake averaging)
# Suitable for mixed effects models with lake as random effect

library(tidyverse)
# Normalize lake names: lowercase, collapse whitespace, remove NBSP and trailing digits
normalize_lake_id <- function(x) {
  x %>%
    stringr::str_replace_all("\u00A0", " ") %>%
    stringr::str_to_lower() %>%
    stringr::str_squish() %>%
    stringr::str_replace_all("[^a-z0-9 ]+", " ") %>%
    stringr::str_squish() %>%
    stringr::str_remove("\\d+$") %>%
    stringr::str_squish() %>%
    stringr::str_replace_all(" ", "")
}

canonicalize_lake_id <- function(x) {
  x <- suppressWarnings(as.character(x))
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    stringr::str_detect(x, "tsiko.*loon|loon.*tsiko|^loon$") ~ "tsiko",
    stringr::str_detect(x, "farewell|^far$") ~ "far",
    stringr::str_detect(x, "rss|rosseau") ~ "rosseau",
    stringr::str_detect(x, "mcnair") ~ "mcn",
    stringr::str_detect(x, "muskeg") ~ "msk",
    stringr::str_detect(x, "sugsaw") ~ "sugsaw",
    stringr::str_detect(x, "grammad|dgramma|^dlake$|^d$") ~ "gramma",
    stringr::str_detect(x, "muchalat") ~ "muchalata",
    stringr::str_detect(x, "ormond") ~ "ormand",
    stringr::str_detect(x, "thiemer|theirmer") ~ "theimer",
    stringr::str_detect(x, "brownbay") ~ "brownbay",
    stringr::str_detect(x, "littlemud") ~ "littlemud",
    stringr::str_detect(x, "littlewoss") ~ "littlewoss",
    stringr::str_detect(x, "lowerstella") ~ "lowerstella",
    stringr::str_detect(x, "uppercampbell") ~ "uppercampbell",
    stringr::str_detect(x, "lowercampbell") ~ "lowercampbell",
    TRUE ~ stringr::str_remove(x, "lake$")
  )
}

# Compute per-lake max depth from sonde and apply 20 m rounding flag rule
compute_sonde_max_depth <- function(sonde_path) {
  sonde <- read.csv(sonde_path, check.names = FALSE)
  if (!"Site" %in% names(sonde) || !"DepthMeters" %in% names(sonde)) {
    stop("Sonde file must contain 'Site' and 'DepthMeters' columns: ", sonde_path)
  }
  sonde %>%
    mutate(
      LakeID = canonicalize_lake_id(normalize_lake_id(Site)),
      DepthMeters_num = suppressWarnings(as.numeric(DepthMeters))
    ) %>%
    group_by(LakeID) %>%
    summarise(
      max_depth_sonde_m = if (all(is.na(DepthMeters_num))) NA_real_ else max(DepthMeters_num, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(flag_sonde_capped_20m = is.finite(max_depth_sonde_m) & round(max_depth_sonde_m, 1) == 20)
}


# Helper: depth-weighted mean via trapezoidal integration
# Falls back to simple mean if <2 valid depth points
depth_weighted_mean <- function(depth_m, values) {
  d <- suppressWarnings(as.numeric(depth_m))
  v <- suppressWarnings(as.numeric(values))
  ok <- !(is.na(d) | is.na(v))
  d <- d[ok]; v <- v[ok]
  if (length(d) == 0 || length(v) == 0) return(NA_real_)
  if (length(d) < 2) return(mean(v, na.rm = TRUE))
  ord <- order(d)
  d <- d[ord]; v <- v[ord]
  total_depth <- max(d) - min(d)
  if (!is.finite(total_depth) || total_depth <= 0) return(mean(v, na.rm = TRUE))
  dz <- diff(d)
  mids <- (v[-1] + v[-length(v)]) / 2
  sum(dz * mids, na.rm = TRUE) / total_depth
}

# Helper: freshwater density (kg/m^3) as a function of temperature (°C)
# Chen and Millero (1977) polynomial approximation for pure water near 1 atm
freshwater_density_kgm3 <- function(temp_C) {
  T <- as.numeric(temp_C)
  rho <- 1000 * (1 - ((T + 288.9414) / (508929.2 * (T + 68.12963))) * (T - 3.9863)^2)
  return(rho)
}

# Compute thermocline depth as depth at maximum buoyancy frequency (maximum stability)
compute_thermocline_depth <- function(depth_m, temp_C) {
  d <- suppressWarnings(as.numeric(depth_m))
  t <- suppressWarnings(as.numeric(temp_C))
  ok <- !(is.na(d) | is.na(t))
  d <- d[ok]; t <- t[ok]
  if (length(d) < 3) return(NA_real_)
  ord <- order(d)
  d <- d[ord]; t <- t[ord]
  rho <- freshwater_density_kgm3(t)
  drho <- diff(rho)
  dz   <- diff(d)
  if (any(!is.finite(dz)) || all(dz == 0)) return(NA_real_)
  g <- 9.81
  rho_mid <- (rho[-1] + rho[-length(rho)]) / 2
  N2 <- -(g / rho_mid) * (drho / dz)  # buoyancy frequency squared (s^-2)
  i_max <- which.max(N2)
  if (length(i_max) == 0 || !is.finite(i_max)) return(NA_real_)
  # Thermocline depth as midpoint of the layer with max N2
  depth_tc <- mean(d[c(i_max, i_max + 1)])
  as.numeric(depth_tc)
}

# Compute depth (m) where DO_mgL crosses threshold (e.g., 1 mg/L) via linear interpolation
compute_do_threshold_depth <- function(depth_m, DO_mgL, threshold = 1) {
  d <- suppressWarnings(as.numeric(depth_m))
  o <- suppressWarnings(as.numeric(DO_mgL))
  ok <- !(is.na(d) | is.na(o))
  d <- d[ok]; o <- o[ok]
  if (length(d) < 2) return(NA_real_)
  ord <- order(d)
  d <- d[ord]; o <- o[ord]
  # Find first segment where it crosses from >thr to <=thr with increasing depth
  for (i in seq_len(length(d) - 1)) {
    o1 <- o[i]; o2 <- o[i + 1]
    if (is.na(o1) || is.na(o2)) next
    if ((o1 > threshold && o2 <= threshold) || (o1 < threshold && o2 >= threshold)) {
      # linear interpolation for depth at threshold
      if (o2 == o1) return(d[i + 1])
      frac <- (threshold - o1) / (o2 - o1)
      z <- d[i] + frac * (d[i + 1] - d[i])
      return(as.numeric(z))
    }
  }
  return(NA_real_)
}

# Compute epilimnion depth (m) as the shallowest depth where |dT/dz| ≥ threshold (°C/m)
# Returns the upper bound of the first segment meeting the threshold (start of metalimnion)
compute_epilimnion_depth <- function(depth_m, temp_C, gradient_threshold = 1) {
  d <- suppressWarnings(as.numeric(depth_m))
  t <- suppressWarnings(as.numeric(temp_C))
  ok <- !(is.na(d) | is.na(t))
  d <- d[ok]; t <- t[ok]
  if (length(d) < 2) return(NA_real_)
  ord <- order(d)
  d <- d[ord]; t <- t[ord]
  dz <- diff(d)
  dt <- diff(t)
  if (any(is.na(dz)) || any(dz <= 0, na.rm = TRUE)) return(NA_real_)
  grad <- abs(dt / dz)
  i <- which(grad >= gradient_threshold)[1]
  if (is.na(i)) return(NA_real_)
  as.numeric(d[i])
}

# Use rLakeAnalyzer to compute thermocline and epilimnion depths after light preprocessing
# - Bin to 0.1 m and average duplicates
# - Remove spike outliers via running median residuals (>2 °C)
compute_layer_depths_rlake <- function(depth_m, temp_C, bin_size_m = 0.1, spike_threshold_C = 2) {
  d <- suppressWarnings(as.numeric(depth_m))
  t <- suppressWarnings(as.numeric(temp_C))
  ok <- !(is.na(d) | is.na(t))
  d <- d[ok]; t <- t[ok]
  if (length(d) < 3) return(list(thermo = NA_real_, epi = NA_real_))
  ord <- order(d)
  d <- d[ord]; t <- t[ord]
  # Bin and average at specified resolution
  d_round <- round(d / bin_size_m) * bin_size_m
  prof <- tibble::tibble(depth = d_round, temp = t) %>%
    dplyr::group_by(depth) %>%
    dplyr::summarise(temp = mean(temp, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(depth)
  if (nrow(prof) < 3) return(list(thermo = NA_real_, epi = NA_real_))
  # Outlier removal via running median residuals
  k <- min(5L, if (nrow(prof) %% 2 == 0) nrow(prof) - 1L else nrow(prof))
  if (is.na(k) || k < 3L) k <- 3L
  med <- stats::runmed(prof$temp, k = k)
  resid <- prof$temp - med
  keep <- is.finite(resid) & abs(resid) <= spike_threshold_C
  prof <- prof[keep, , drop = FALSE]
  if (nrow(prof) < 3) return(list(thermo = NA_real_, epi = NA_real_))
  # rLakeAnalyzer if available
  if (requireNamespace("rLakeAnalyzer", quietly = TRUE)) {
    thermo <- suppressWarnings(try(rLakeAnalyzer::thermo.depth(wtr = prof$temp, depths = prof$depth), silent = TRUE))
    meta <- suppressWarnings(try(rLakeAnalyzer::meta.depths(wtr = prof$temp, depths = prof$depth), silent = TRUE))
    thermo_val <- if (inherits(thermo, "try-error")) NA_real_ else as.numeric(thermo)
    epi_val <- NA_real_
    if (!inherits(meta, "try-error")) {
      # meta.depths returns c(top, bottom)
      epi_val <- suppressWarnings(as.numeric(meta[1]))
    }
    return(list(thermo = thermo_val, epi = epi_val))
  }
  # If rLakeAnalyzer is not available, return NAs so issues are visible
  return(list(thermo = NA_real_, epi = NA_real_))
}

# Function 1: Process Sonde Data - Group by Deployment and Classify as Shallow/Deep
process_sonde_deployments <- function(sonde_path) {
  message("Processing sonde data by deployment groups...")
  
  # Read data
  sonde_df <- read.csv(sonde_path)
  
  # Debug: Check original Site column
  message("Original Site examples:")
  original_sites <- sonde_df %>%
    dplyr::select(Site) %>%
    distinct() %>%
    head(10)
  print(original_sites)
  
  # Separate site into lake ID and sample number
  data_separated <- sonde_df %>%
    separate(Site, into = c("LakeID", "SampleNumber"), 
             sep = "(?<=[a-zA-Z])(?=\\d)",
             fill = "right", remove = FALSE)  # Keep original Site column
  
  # Debug: Check site separation
  message("Site separation examples:")
  site_check <- data_separated %>%
    dplyr::select(Site, LakeID, SampleNumber) %>%
    distinct() %>%
    head(10)
  print(site_check)
  
  # Handle missing sample numbers and standardize lake names
  sonde_clean <- data_separated %>%
    mutate(SampleNumber = if_else(is.na(SampleNumber), "1", SampleNumber)) %>%
    mutate(LakeID = canonicalize_lake_id(normalize_lake_id(LakeID))) %>%
    filter(DepthMeters >= 0) %>%
    # Remove unwanted columns
    dplyr::select(-any_of(c("User.ID", "CableVNotSure"))) %>%
    # Rename for consistency
    mutate(
      Depth_m = DepthMeters,
      Temp = TemperatureCelcius,
      DO_percent = PercentDissolvedOxygen,
      DO_mgL = DissolvedOxygen_mg_L,
      Phyco = PhycocyaninRFU,
      Chl = ChlorophyllRFU,
      pH = pH,
      fDOM = fDOMRFU
    ) %>%
    # Add time parsing if available
    mutate(
      Time_parsed = if("Time" %in% names(.)) {
        as.POSIXct(paste(Date, Time), format="%d-%b-%y %H:%M", tz="UTC")
      } else {
        as.POSIXct(NA)
      }
    ) %>%
    arrange(LakeID, Time_parsed, Depth_m)
  
  # Group sonde measurements into deployments
  # Step 1: Check which lakes have explicit replicates
  lakes_with_replicates <- sonde_clean %>%
    group_by(LakeID) %>%
    summarise(has_explicit_replicates = length(unique(SampleNumber)) > 1, .groups = "drop")
  
  # Step 2: Process each lake separately based on replicate status
  sonde_deployments <- sonde_clean %>%
    left_join(lakes_with_replicates, by = "LakeID") %>%
    group_by(LakeID) %>%
    arrange(Time_parsed, Depth_m) %>%
    mutate(
      # For lakes WITHOUT explicit replicates, detect deployment breaks
      # More sensitive depth reset: if depth decreases significantly AND we were deep (>5m)
      was_deep = lag(Depth_m, default = 0) > 5,
      depth_decrease = c(0, diff(Depth_m)),
      depth_reset = ifelse(!has_explicit_replicates, 
                          was_deep & (depth_decrease < -3), FALSE),  # More sensitive
      # Time gap (reduced threshold to catch 4-minute gaps)
      time_gap = ifelse(!has_explicit_replicates & !all(is.na(Time_parsed)), 
                       c(FALSE, as.numeric(diff(Time_parsed, units="mins")) > 3), FALSE),
      # Combined break detection
      deployment_break = ifelse(!has_explicit_replicates, 
                               depth_reset | time_gap | (row_number() == 1), FALSE),
      # Create sub-deployment groups
      auto_deployment = cumsum(deployment_break),
      # Final deployment ID
      deployment_id = ifelse(has_explicit_replicates,
                            paste(LakeID, SampleNumber, sep = "_"),
                            paste(LakeID, auto_deployment, sep = "_"))
    ) %>%
    ungroup() %>%
    # Group by final deployment ID for summary stats
    group_by(deployment_id) %>%
    mutate(
      max_depth_deployment = max(Depth_m, na.rm = TRUE),
      min_depth_deployment = min(Depth_m, na.rm = TRUE),
      n_measurements = n()
    ) %>%
    ungroup()
  
  # Debug: Check how many deployments per lake
  deployments_per_lake <- sonde_deployments %>%
    dplyr::select(LakeID, SampleNumber, deployment_id) %>%
    distinct() %>%
    count(LakeID, name = "n_deployments")
  
  message("Deployments per lake (first 10):")
  print(head(deployments_per_lake, 10))
  
  # Special debug for cedar
  cedar_debug <- sonde_deployments %>%
    filter(LakeID == "cedar") %>%
    dplyr::select(LakeID, Depth_m, Time_parsed, was_deep, depth_decrease, depth_reset, 
                  time_gap, deployment_break, deployment_id) %>%
    head(20)
  
  if(nrow(cedar_debug) > 0) {
    message("\nCedar deployment detection debug:")
    print(cedar_debug)
  }
  
  # Classify deployments as shallow or deep
  deployment_summary <- sonde_deployments %>%
    group_by(LakeID, deployment_id) %>%
    summarise(
      max_depth = max(Depth_m, na.rm = TRUE),
      min_depth = min(Depth_m, na.rm = TRUE),
      n_points = n(),
      .groups = "drop"
    ) %>%
    group_by(LakeID) %>%
    mutate(
      # Simple relative classification: deepest deployment = "deep", others = "shallow"
      deployment_type = case_when(
        n() == 1 ~ "single",
        max_depth == max(max_depth) ~ "deep",    # The deepest deployment
        TRUE ~ "shallow"                         # All other deployments
      )
    ) %>%
    ungroup()
  
  # Calculate depth-integrated (depth-weighted) averages for each deployment
  deployment_conditions <- sonde_deployments %>%
    left_join(deployment_summary %>% dplyr::select(LakeID, deployment_id, deployment_type, max_depth), 
              by = c("LakeID", "deployment_id")) %>%
    group_by(LakeID, deployment_id, deployment_type) %>%
    summarise(
      # Depth-weighted (trapezoidal) conditions over full profile
      integrated_temp = depth_weighted_mean(Depth_m, Temp),
      integrated_pH = depth_weighted_mean(Depth_m, pH),
      integrated_DO_percent = depth_weighted_mean(Depth_m, pmax(0, DO_percent)),
      integrated_DO_mgL = depth_weighted_mean(Depth_m, pmax(0, DO_mgL)),
      integrated_chl = depth_weighted_mean(Depth_m, Chl),
      integrated_phyco = depth_weighted_mean(Depth_m, Phyco),
      integrated_fDOM = depth_weighted_mean(Depth_m, fDOM),

      # Thermocline, oxygen threshold (2 mg/L), and epilimnion depths
      thermocline_depth_m = compute_layer_depths_rlake(Depth_m, Temp)$thermo,
      oxygen_threshold_depth_m = compute_do_threshold_depth(Depth_m, pmax(0, DO_mgL), threshold = 2),
      epilimnion_depth_m = compute_layer_depths_rlake(Depth_m, Temp)$epi,
      refuge_bottom_depth_m = ifelse(is.na(oxygen_threshold_depth_m), max(Depth_m, na.rm = TRUE), oxygen_threshold_depth_m),
      refuge_is_bottom_bounded = is.na(oxygen_threshold_depth_m),
      refuge_thickness_m = {
        top <- epilimnion_depth_m
        bot <- refuge_bottom_depth_m
        if (is.na(top) || is.na(bot)) NA_real_ else max(0, bot - top)
      },
      refuge_m = refuge_thickness_m,
      
      # Surface conditions (top 2m), depth-weighted within window when possible
      surface_temp = depth_weighted_mean(Depth_m[Depth_m <= 2], Temp[Depth_m <= 2]),
      surface_pH = depth_weighted_mean(Depth_m[Depth_m <= 2], pH[Depth_m <= 2]),
      surface_DO_percent = depth_weighted_mean(Depth_m[Depth_m <= 2], pmax(0, DO_percent[Depth_m <= 2])),
      surface_chl = depth_weighted_mean(Depth_m[Depth_m <= 2], Chl[Depth_m <= 2]),
      
      # Bottom conditions (deepest 20% of profile), depth-weighted
      bottom_temp = {
        cut <- suppressWarnings(quantile(Depth_m, 0.8, na.rm = TRUE))
        depth_weighted_mean(Depth_m[Depth_m >= cut], Temp[Depth_m >= cut])
      },
      bottom_DO_percent = {
        cut <- suppressWarnings(quantile(Depth_m, 0.8, na.rm = TRUE))
        depth_weighted_mean(Depth_m[Depth_m >= cut], pmax(0, DO_percent[Depth_m >= cut]))
      },
      bottom_DO_mgL = {
        cut <- suppressWarnings(quantile(Depth_m, 0.8, na.rm = TRUE))
        depth_weighted_mean(Depth_m[Depth_m >= cut], pmax(0, DO_mgL[Depth_m >= cut]))
      },
      min_DO_mgL = suppressWarnings(pmax(0, min(DO_mgL, na.rm = TRUE))),
      
      # Deployment metadata
      max_depth_sampled = max(Depth_m, na.rm = TRUE),
      n_depth_points = n(),
      depth_range = max(Depth_m, na.rm = TRUE) - min(Depth_m, na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    # Create deployment-specific ID for matching
    mutate(DeploymentID = paste(LakeID, deployment_type, sep = "_"))
  
  message("Sonde deployments processed:")
  message("- ", nrow(deployment_conditions), " deployments identified")
  message("- ", sum(deployment_conditions$deployment_type == "shallow"), " shallow deployments")
  message("- ", sum(deployment_conditions$deployment_type == "deep"), " deep deployments")
  message("- ", sum(deployment_conditions$deployment_type == "single"), " single deployments")
  
  # Debug: show a few deployment examples
  message("\nExample deployments:")
  example_deployments <- deployment_conditions %>% 
    dplyr::select(LakeID, deployment_type, max_depth_sampled, n_depth_points, integrated_temp, integrated_DO_percent) %>%
    head(5)
  print(example_deployments)
  
  # Return both detailed data and deployment summaries
  result <- list(
    deployment_conditions = deployment_conditions,
    raw_deployments = sonde_deployments,
    deployment_summary = deployment_summary
  )
  
  return(result)
}

# Function 2: Process Zooplankton Data at Sample Level
process_zooplankton_samples <- function(zoop_path) {
  message("Processing zooplankton data at sample level...")
  
  # Read data
  zoop_df <- read.csv(zoop_path)
  
  # Parameters
  subsample_volume_ml <- 2
  net_radius_m <- 0.127 / 2
  NumberOfTows <- 2
  
  # Define taxon columns
  taxon_cols <- c("Calanoids", "Cyclopoids", "Nauplii", "Daphnia", "Bosmina", 
                  "Chydorus", "Diaphanosoma", "Holopedium", "Chaoborus", 
                  "WaterMites", "Brachionus.Keratella", "ConochilusColonies", 
                  "Asplanchna", "Ciliate", "Amoeboid", "Rotifers..other.", 
                  "Polyphemidae")
  
  # Keep only existing columns
  taxon_cols <- taxon_cols[taxon_cols %in% names(zoop_df)]
  
  # Convert to numeric
  zoop_df[taxon_cols] <- lapply(zoop_df[taxon_cols], function(x) as.numeric(as.character(x)))
  
  # Process each sample
  zoop_samples <- zoop_df %>%
    rename(LakeID = LakeName) %>%
    mutate(
      Depth_Tow_m = as.numeric(Depth_Tow_m),
      # Normalize LakeID to align with sonde/geography naming
      LakeID = canonicalize_lake_id(normalize_lake_id(LakeID)),
      # Create unique sample ID
      SampleID = paste(LakeID, row_number(), sep = "_zoop_")
    ) %>%
    filter(!is.na(Depth_Tow_m))
  
  # Calculate densities for each sample in long format
  zoop_long <- zoop_samples %>%
    pivot_longer(cols = all_of(taxon_cols), names_to = "Taxon", values_to = "Count") %>%
    mutate(
      Count = as.numeric(Count),
      counted_volume_ml = No2mlSubsamplesCounted * subsample_volume_ml,
      total_individuals_in_tube = (Count * DilutionVolume) / counted_volume_ml,
      volume_sampled_m3 = pi * net_radius_m^2 * Depth_Tow_m,
      volume_sampled_L = volume_sampled_m3 * 1000,
      density_per_m3 = (total_individuals_in_tube / volume_sampled_m3) / NumberOfTows,
      density_per_L = (total_individuals_in_tube / volume_sampled_L) / NumberOfTows
    ) %>%
    filter(!is.na(Count))
  
  # Pivot wider to get one row per sample
  zoop_wide <- zoop_long %>%
    dplyr::select(SampleID, LakeID, Depth_Tow_m, Taxon, density_per_L) %>%
    pivot_wider(
      names_from = Taxon,
      values_from = density_per_L,
      names_glue = "{Taxon}_{.value}"
    )
  
  # Add sample metadata
  # Check which columns exist
  metadata_cols <- c("SampleID", "LakeID", "Depth_Tow_m", "DilutionVolume", 
                     "No2mlSubsamplesCounted")
  optional_cols <- c("Date", "Time")
  
  # Add optional columns if they exist
  available_cols <- c(metadata_cols, optional_cols[optional_cols %in% names(zoop_samples)])
  
  # Ensure LakeID exists in both datasets
  if(!"LakeID" %in% names(zoop_samples)) {
    stop("LakeID column missing from zooplankton data")
  }
  if(!"LakeID" %in% names(zoop_wide)) {
    stop("LakeID column missing from zooplankton wide data")
  }
  
  sample_metadata <- zoop_samples %>%
    dplyr::select(all_of(available_cols)) %>%
    distinct()
  
  # Check join keys exist
  join_keys <- c("SampleID", "LakeID", "Depth_Tow_m")
  join_keys <- join_keys[join_keys %in% names(sample_metadata) & join_keys %in% names(zoop_wide)]
  
  zoop_final <- sample_metadata %>%
    left_join(zoop_wide, by = join_keys)
  
  # === Add per-sample diversity metrics (exclude Nauplii) ===
  density_cols <- names(zoop_final)[grepl("_density_per_L$", names(zoop_final))]
  density_cols <- setdiff(density_cols, c("Nauplii_density_per_L"))
  if (length(density_cols) > 0) {
    dens_mat <- as.matrix(zoop_final[ , density_cols, drop = FALSE])
    dens_mat[is.na(dens_mat)] <- 0
    row_totals <- rowSums(dens_mat)
    with_positive_total <- row_totals > 0 & is.finite(row_totals)
    p_mat <- matrix(0, nrow = nrow(dens_mat), ncol = ncol(dens_mat))
    if (any(with_positive_total)) {
      p_mat[with_positive_total, ] <- dens_mat[with_positive_total, , drop = FALSE] / row_totals[with_positive_total]
    }
    sh_terms <- p_mat
    sh_terms[sh_terms <= 0] <- NA_real_
    diversity_shannon <- -rowSums(sh_terms * log(sh_terms), na.rm = TRUE)
    D <- rowSums(p_mat^2, na.rm = TRUE)
    diversity_simpson <- 1 - D
    diversity_invsimpson <- ifelse(D > 0, 1 / D, NA_real_)
    diversity_richness <- rowSums(dens_mat > 0, na.rm = TRUE)
    diversity_shannon[!with_positive_total] <- NA_real_
    diversity_simpson[!with_positive_total] <- NA_real_
    diversity_invsimpson[!with_positive_total] <- NA_real_
    zoop_final <- zoop_final %>% mutate(
      diversity_shannon = diversity_shannon,
      diversity_simpson = diversity_simpson,
      diversity_invsimpson = diversity_invsimpson,
      diversity_richness = diversity_richness
    )
  }
  
  return(zoop_final)
}

# Function 2b: Process Copepod Infection Data 
# Expects columns like: LAKE_CODE, LAKE, DEPTH (DEEP/SHALLOW), MONTH, YEAR, ... infection metrics
process_infection_data <- function(infection_path) {
  message("Processing infection data...")

  inf_raw <- read.csv(infection_path, check.names = FALSE)
  # Trim/normalize column names (the file has e.g. ' YEAR' with a leading space)
  names(inf_raw) <- stringr::str_squish(names(inf_raw))

  required <- c("LAKE", "DEPTH")
  missing_req <- setdiff(required, names(inf_raw))
  if (length(missing_req) > 0) {
    stop("Infection file missing required columns: ", paste(missing_req, collapse = ", "), " in ", infection_path)
  }

  inf <- inf_raw %>%
    mutate(
      LAKE = as.character(LAKE),
      LAKE_CODE = if ("LAKE_CODE" %in% names(.)) as.character(LAKE_CODE) else NA_character_,
      MONTH = if ("MONTH" %in% names(.)) as.character(MONTH) else NA_character_,
      YEAR = if ("YEAR" %in% names(.)) suppressWarnings(as.integer(YEAR)) else NA_integer_,
      DEPTH = as.character(DEPTH)
    ) %>%
    mutate(
      LakeID = canonicalize_lake_id(normalize_lake_id(LAKE)),
      infection_sample_type = case_when(
        tolower(trimws(DEPTH)) %in% c("deep", "d") ~ "deep",
        tolower(trimws(DEPTH)) %in% c("shallow", "s") ~ "shallow",
        TRUE ~ NA_character_
      ),
      MatchID = paste(LakeID, infection_sample_type, sep = "_")
    )

  # Keep only infection metrics that exist; prefix to avoid clobbering main columns
  col_map <- c(
    "RAW_COPE_COPIES" = "infection_raw_cope_copies",
    "ADJUSTED_COPE_COPIES" = "infection_adjusted_cope_copies",
    "COPE_COUNTS" = "infection_cope_counts",
    "RAW_WORM_COPIES" = "infection_raw_worm_copies",
    "ADJUSTED_WORM_COPIES" = "infection_adjusted_worm_copies",
    "WORM_COUNTS" = "infection_worm_counts",
    "INFECTION_BURDEN" = "infection_burden"
  )

  present <- names(col_map)[names(col_map) %in% names(inf)]
  inf_out <- inf %>%
    transmute(
      MatchID,
      LakeID,
      infection_sample_type,
      infection_lake_code = LAKE_CODE,
      infection_month = MONTH,
      infection_year = YEAR,
      !!!rlang::set_names(lapply(present, function(nm) suppressWarnings(as.numeric(inf[[nm]]))), col_map[present])
    ) %>%
    distinct()

  # Guard against duplicated keys (would explode row count on join)
  dup_keys <- inf_out %>% count(MatchID) %>% filter(n > 1)
  if (nrow(dup_keys) > 0) {
    stop(
      "Infection data has duplicated (LakeID, DEPTH) keys after normalization. Duplicates:\n",
      paste(dup_keys$MatchID, collapse = ", ")
    )
  }

  return(inf_out)
}

# Function 3: Get Lake-Level Variables (constant within lake)
get_lake_level_data <- function(bathy_path, geo_path, sonde_path = "data/Sonde Data/VancouverSondeData2023_UTF8.csv") {
  message("Processing lake-level variables...")
  
  # Process bathymetry - support CSV or Excel based on extension
  bathy_df <- NULL
  if (grepl("\\.xlsx?$", bathy_path, ignore.case = TRUE)) {
    if(!requireNamespace("readxl", quietly = TRUE)) {
      stop("readxl package is required to read Excel. install.packages('readxl')")
    }
    bathy_df <- readxl::read_excel(bathy_path)
  } else if (grepl("\\.csv$", bathy_path, ignore.case = TRUE)) {
    bathy_df <- read.csv(bathy_path, check.names = FALSE)
  } else {
    stop("Unsupported bathymetry format: ", bathy_path)
  }
  
  if("LakeID" %in% names(bathy_df)) {
    lake_bathy <- bathy_df
  } else if("Lake" %in% names(bathy_df)) {
    lake_bathy <- bathy_df %>% rename(LakeID = Lake)
  } else {
    stop("No Lake or LakeID column found in bathymetry data")
  }
  
  # Check for required columns and rename if needed
  col_mapping <- c(
    "MaxDepth_m" = "MaxDepth",
    "Surface_Area_ha" = "SurfaceArea_ha",
    "MeanDepth_m" = "MeanDepth"
  )
  
  for(old_name in names(col_mapping)) {
    if(old_name %in% names(lake_bathy) && !col_mapping[old_name] %in% names(lake_bathy)) {
      lake_bathy <- lake_bathy %>% rename(!!col_mapping[old_name] := !!old_name)
    }
  }
  
  # Process numeric columns
  numeric_cols <- c("MaxDepth", "MeanDepth", "SurfaceArea_ha", "Depth_Ratio")
  for(col in numeric_cols) {
    if(col %in% names(lake_bathy)) {
      lake_bathy[[col]] <- as.numeric(gsub("n/a", NA, lake_bathy[[col]]))
    }
  }
  # Normalize LakeID early
  lake_bathy <- lake_bathy %>% mutate(LakeID = canonicalize_lake_id(normalize_lake_id(LakeID)))

  # === Fill MaxDepth using sonde when bathymetry missing; flag 20 m rounded sonde as missing ===
  sonde_max <- try(compute_sonde_max_depth(sonde_path), silent = TRUE)
  if (inherits(sonde_max, "try-error")) {
    warning("Failed to compute sonde max depth: ", as.character(sonde_max))
    sonde_max <- tibble::tibble(LakeID = character(), max_depth_sonde_m = numeric(), flag_sonde_capped_20m = logical())
  }

  lake_bathy <- lake_bathy %>%
    left_join(sonde_max, by = "LakeID") %>%
    mutate(
      MaxDepth_bathy = MaxDepth,
      MaxDepth_source = dplyr::case_when(
        is.finite(MaxDepth_bathy) ~ "bathymetry",
        is.finite(max_depth_sonde_m) & !flag_sonde_capped_20m ~ "sonde",
        TRUE ~ "missing"
      ),
      MaxDepth = dplyr::case_when(
        is.finite(MaxDepth_bathy) ~ MaxDepth_bathy,
        is.finite(max_depth_sonde_m) & !flag_sonde_capped_20m ~ max_depth_sonde_m,
        TRUE ~ NA_real_
      ),
      MaxDepth_lower_bound_m = dplyr::if_else(!is.finite(MaxDepth) & flag_sonde_capped_20m, 20, NA_real_),
      MaxDepth_missing_both = is.na(MaxDepth_bathy) & is.na(max_depth_sonde_m)
    )
  
  # Calculate derived metrics - check which columns exist
  lake_bathy <- lake_bathy %>%
    mutate(
      # Only calculate if column exists
      SurfaceArea_m2 = if("SurfaceArea_ha" %in% names(.)) SurfaceArea_ha * 10000 else NA_real_,
      Volume_m3 = if(all(c("SurfaceArea_m2", "MeanDepth") %in% names(.))) SurfaceArea_m2 * MeanDepth else NA_real_,
      depth_ratio = case_when(
        "Depth_Ratio" %in% names(.) & !is.na(Depth_Ratio) ~ Depth_Ratio,
        all(c("MeanDepth", "MaxDepth") %in% names(.)) ~ MeanDepth / MaxDepth,
        TRUE ~ NA_real_
      )
    )
  # Invalidate impossible or degenerate values: depth_ratio outside (0,1)
  lake_bathy <- lake_bathy %>%
    mutate(
      depth_ratio = ifelse(!is.finite(depth_ratio) | depth_ratio <= 0 | depth_ratio >= 1, NA_real_, depth_ratio),
      # If MeanDepth > MaxDepth, set MeanDepth and derived Volume_m3 to NA so downstream doesn't use it
      MeanDepth = ifelse(is.finite(MeanDepth) & is.finite(MaxDepth) & MeanDepth > MaxDepth, NA_real_, MeanDepth)
    )
  # Ensure volume is computed when inputs available
  lake_bathy <- lake_bathy %>%
    mutate(
      # Impute MeanDepth when depth_ratio and MaxDepth exist
      MeanDepth = ifelse(is.na(MeanDepth) & is.finite(depth_ratio) & is.finite(MaxDepth), depth_ratio * MaxDepth, MeanDepth),
      Volume_m3 = ifelse(!is.na(SurfaceArea_m2) & !is.na(MeanDepth), SurfaceArea_m2 * MeanDepth, Volume_m3)
    )
  
  # Select only columns that exist (include provenance and flags)
  available_cols <- c(
    "LakeID", "MaxDepth", "SurfaceArea_m2", "Volume_m3", "MeanDepth", "depth_ratio",
    "MaxDepth_bathy", "max_depth_sonde_m", "MaxDepth_source", "MaxDepth_lower_bound_m", "MaxDepth_missing_both"
  )
  available_cols <- available_cols[available_cols %in% names(lake_bathy)]
  
  lake_bathy <- lake_bathy %>%
    dplyr::select(all_of(available_cols))
  
  # Process geography
  geo_df <- read.csv(geo_path)
  
  lake_geo <- geo_df %>%
    rename(LakeID = names(.)[1]) %>%
    mutate(
      distance_from_ocean_m = as.numeric(distance.from.ocean..m.),
      elevation_m = as.numeric(elevation..m.),
      longitude = as.numeric(longitude),
      latitude = as.numeric(latitude)
    ) %>%
    mutate(
      LakeID_clean = normalize_lake_id(LakeID),
      LakeID = case_when(
        grepl("nowgs|experiment|^\\*", LakeID_clean, ignore.case = TRUE) ~ NA_character_,
        LakeID_clean == "" ~ NA_character_,
        TRUE ~ canonicalize_lake_id(LakeID_clean)
      )
    ) %>%
    filter(!is.na(LakeID)) %>%
    dplyr::select(
      LakeID,
      Ocean_Distance_m = distance_from_ocean_m,
      Lake_Level_m = elevation_m,
      Latitude = latitude,
      Longitude = longitude
    )
  
  # Join lake-level data
  lake_level_vars <- lake_bathy %>%
    left_join(lake_geo, by = "LakeID") %>%
    rename(DFO_m = Ocean_Distance_m)

  # === Optional: Join Land Use data if available ===
  landuse_path <- "data/LandUse/land_use.csv"
  if (file.exists(landuse_path)) {
    land_use <- read.csv(landuse_path, check.names = FALSE) %>%
      # Expect a column named 'lake_name' with values like "X Lake"
      dplyr::rename(lake_name = lake_name) %>%
      dplyr::mutate(
        LakeID_clean = normalize_lake_id(lake_name),
        LakeID = dplyr::case_when(
          grepl("nowgs|experiment|^\\*", LakeID_clean, ignore.case = TRUE) ~ NA_character_,
          LakeID_clean == "" ~ NA_character_,
          TRUE ~ canonicalize_lake_id(LakeID_clean)
        )
      ) %>%
      dplyr::filter(!is.na(LakeID)) %>%
      dplyr::transmute(
        LakeID,
        landuse_avg_chlorophyll = suppressWarnings(as.numeric(avg_chlorophyll)),
        landuse_mean_cover_pattern = suppressWarnings(as.numeric(mean_cover_pattern)),
        landuse_mean_crown_closure = suppressWarnings(as.numeric(mean_crown_closure)),
        landuse_elevation = suppressWarnings(as.numeric(elevation)),
        landuse_mean_height = suppressWarnings(as.numeric(mean_height)),
        landuse_loose_road_density = suppressWarnings(as.numeric(loose_road_density)),
        landuse_paved_road_density = suppressWarnings(as.numeric(paved_road_density)),
        landuse_mean_forest_age = suppressWarnings(as.numeric(mean_forest_age))
      )

    lake_level_vars <- lake_level_vars %>%
      dplyr::left_join(land_use, by = "LakeID")
  }
  
  return(lake_level_vars)
}

# Main function to create sample-level dataset with deployment-matched sonde data
create_sample_level_dataframe <- function(
  sonde_path = "data/Sonde Data/VancouverSondeData2023_UTF8.csv",
  bathy_path = "data/Bathymetry Data/LakeBathymetry.csv",
  geo_path = "data/LakeGeography/Field_2023_Lake_Level_Ocean_Distance.csv",
  zoop_path = "data/Zooplankton ID Data/ZoopIDDataUpdated.csv",
  infection_path = "data/Infection Data/42_Lake_Copepod_Infection.csv"
) {
  
  # Process each dataset
  sonde_result <- process_sonde_deployments(sonde_path)
  zoop_samples <- process_zooplankton_samples(zoop_path)
  lake_vars <- get_lake_level_data(bathy_path, geo_path)
  
  # Debug: Check column names
  message("zoop_samples columns: ", paste(names(zoop_samples)[1:10], collapse = ", "))
  message("lake_vars columns: ", paste(names(lake_vars)[1:5], collapse = ", "))
  
  # Join zooplankton and lake-level data
  # Check if LakeID exists in both datasets
  if(!"LakeID" %in% names(zoop_samples)) {
    stop("LakeID column missing from zoop_samples. Available: ", paste(names(zoop_samples), collapse = ", "))
  }
  if(!"LakeID" %in% names(lake_vars)) {
    stop("LakeID column missing from lake_vars. Available: ", paste(names(lake_vars), collapse = ", "))
  }
  
  zoop_with_lake <- zoop_samples %>%
    left_join(lake_vars, by = "LakeID")
  
  message("\nMatching sonde deployments to zooplankton samples...")
  
  # Classify zooplankton samples as shallow or deep
  zoop_classified <- zoop_with_lake %>%
    group_by(LakeID) %>%
    mutate(
      n_samples_lake = n(),
      depth_rank = rank(-Depth_Tow_m, na.last = "keep"),
      zoop_sample_type = case_when(
        n_samples_lake == 1 ~ "single",
        depth_rank == 1 ~ "deep",
        TRUE ~ "shallow"
      )
    ) %>%
    ungroup() %>%
    # Create matching key
    mutate(
      MatchID = paste(LakeID, zoop_sample_type, sep = "_"),
      MatchID_single = paste(LakeID, "single", sep = "_")
    )
  
  # Debug: Check zoop_classified columns
  message("zoop_classified columns: ", paste(names(zoop_classified)[1:10], collapse = ", "))
  message("deployment_conditions columns: ", paste(names(sonde_result$deployment_conditions)[1:10], collapse = ", "))
  
  # Match zooplankton to sonde deployments
  matched_data <- zoop_classified %>%
    left_join(
      sonde_result$deployment_conditions %>%
        rename(
          MatchID = DeploymentID,
          sonde_deployment_type = deployment_type
        ) %>%
        dplyr::select(-LakeID),  # Remove LakeID to avoid .x/.y issue
      by = "MatchID"
    )
  # If no exact deep/shallow match, coalesce from the lake's single deployment (if present)
  single_dc <- sonde_result$deployment_conditions %>%
    filter(deployment_type == "single") %>%
    mutate(MatchID_single = paste(LakeID, "single", sep = "_")) %>%
    dplyr::select(-LakeID, -DeploymentID) %>%
    rename(deployment_type_single = deployment_type)

  matched_data <- matched_data %>%
    left_join(single_dc, by = "MatchID_single", suffix = c("", "_single")) %>%
    mutate(
      thermocline_depth_m = dplyr::coalesce(thermocline_depth_m, thermocline_depth_m_single),
      epilimnion_depth_m = dplyr::coalesce(epilimnion_depth_m, epilimnion_depth_m_single),
      oxygen_threshold_depth_m = dplyr::coalesce(oxygen_threshold_depth_m, oxygen_threshold_depth_m_single),
      integrated_temp = dplyr::coalesce(integrated_temp, integrated_temp_single),
      integrated_pH = dplyr::coalesce(integrated_pH, integrated_pH_single),
      integrated_DO_percent = dplyr::coalesce(integrated_DO_percent, integrated_DO_percent_single),
      integrated_DO_mgL = dplyr::coalesce(integrated_DO_mgL, integrated_DO_mgL_single),
      integrated_chl = dplyr::coalesce(integrated_chl, integrated_chl_single),
      integrated_phyco = dplyr::coalesce(integrated_phyco, integrated_phyco_single),
      integrated_fDOM = dplyr::coalesce(integrated_fDOM, integrated_fDOM_single),
      sonde_deployment_type = dplyr::coalesce(sonde_deployment_type, ifelse(!is.na(deployment_type_single), "single", NA_character_))
    ) %>%
    dplyr::select(-dplyr::ends_with("_single"), -deployment_type_single)

  # Fallback hierarchy for structure metrics: exact match -> single -> deep -> shallow
  # Build per-lake deep and shallow summaries (first non-NA per lake, no numeric averaging)
  pick_first <- function(x) {
    ix <- which(!is.na(x))[1]
    if (length(ix) == 0 || is.na(ix)) return(NA_real_)
    x[ix]
  }
  deep_dc <- sonde_result$deployment_conditions %>%
    filter(deployment_type == "deep") %>%
    group_by(LakeID) %>%
    summarise(
      thermocline_depth_m_deep = pick_first(thermocline_depth_m),
      epilimnion_depth_m_deep = pick_first(epilimnion_depth_m),
      oxygen_threshold_depth_m_deep = pick_first(oxygen_threshold_depth_m),
      integrated_temp_deep = pick_first(integrated_temp),
      integrated_pH_deep = pick_first(integrated_pH),
      integrated_DO_percent_deep = pick_first(integrated_DO_percent),
      integrated_DO_mgL_deep = pick_first(integrated_DO_mgL),
      integrated_chl_deep = pick_first(integrated_chl),
      integrated_phyco_deep = pick_first(integrated_phyco),
      integrated_fDOM_deep = pick_first(integrated_fDOM),
      .groups = "drop"
    )
  matched_data <- matched_data %>%
    left_join(deep_dc, by = "LakeID") %>%
    mutate(
      thermocline_depth_m = dplyr::coalesce(thermocline_depth_m, thermocline_depth_m_deep),
      epilimnion_depth_m = dplyr::coalesce(epilimnion_depth_m, epilimnion_depth_m_deep),
      oxygen_threshold_depth_m = dplyr::coalesce(oxygen_threshold_depth_m, oxygen_threshold_depth_m_deep),
      integrated_temp = dplyr::coalesce(integrated_temp, integrated_temp_deep),
      integrated_pH = dplyr::coalesce(integrated_pH, integrated_pH_deep),
      integrated_DO_percent = dplyr::coalesce(integrated_DO_percent, integrated_DO_percent_deep),
      integrated_DO_mgL = dplyr::coalesce(integrated_DO_mgL, integrated_DO_mgL_deep),
      integrated_chl = dplyr::coalesce(integrated_chl, integrated_chl_deep),
      integrated_phyco = dplyr::coalesce(integrated_phyco, integrated_phyco_deep),
      integrated_fDOM = dplyr::coalesce(integrated_fDOM, integrated_fDOM_deep),
      sonde_deployment_type = dplyr::coalesce(sonde_deployment_type,
        ifelse(!is.na(thermocline_depth_m_deep) | !is.na(integrated_temp_deep), "deep", NA_character_))
    ) %>%
    dplyr::select(-dplyr::matches("_deep$"))

  
  message("After first join, columns: ", paste(names(matched_data)[1:15], collapse = ", "))
  
  # Create fallback data
  fallback_data <- sonde_result$deployment_conditions %>%
    group_by(LakeID) %>%
    summarise(
      fallback_temp = mean(integrated_temp, na.rm = TRUE),
      fallback_pH = mean(integrated_pH, na.rm = TRUE),
      fallback_DO_percent = mean(integrated_DO_percent, na.rm = TRUE),
      fallback_chl = mean(integrated_chl, na.rm = TRUE),
      .groups = "drop"
    )
  
  message("fallback_data columns: ", paste(names(fallback_data), collapse = ", "))
  
  # Add fall-back lake averages for unmatched samples
  matched_data <- matched_data %>%
    left_join(fallback_data, by = "LakeID") %>%
    # Use deployment-specific values when available, otherwise use lake averages
    mutate(
      final_temp = if_else(is.na(integrated_temp), fallback_temp, integrated_temp),
      final_pH = if_else(is.na(integrated_pH), fallback_pH, integrated_pH),
      final_DO_percent = if_else(is.na(integrated_DO_percent), fallback_DO_percent, integrated_DO_percent),
      final_chl = if_else(is.na(integrated_chl), fallback_chl, integrated_chl),
      
      # Quality flags
      sonde_match_quality = case_when(
        !is.na(integrated_temp) ~ "deployment_matched",
        !is.na(fallback_temp) ~ "lake_average",
        TRUE ~ "no_sonde_data"
      )
    )
  
  message("Sonde-zooplankton matching completed:")
  message("- ", sum(matched_data$sonde_match_quality == "deployment_matched", na.rm=TRUE), 
          " samples matched to specific deployments")
  message("- ", sum(matched_data$sonde_match_quality == "lake_average", na.rm=TRUE), 
          " samples using lake averages")
  message("- ", sum(matched_data$sonde_match_quality == "no_sonde_data", na.rm=TRUE), 
          " samples with no sonde data")
  
  # Debug: Check which samples are using lake averages
  lake_avg_samples <- matched_data %>%
    filter(sonde_match_quality == "lake_average") %>%
    dplyr::select(LakeID, zoop_sample_type, MatchID) %>%
    distinct()
  
  message("\nSamples using lake averages:")
  print(lake_avg_samples)
  
  # Debug: Check what sonde deployments are available
  available_deployments <- sonde_result$deployment_conditions %>%
    dplyr::select(LakeID, deployment_type, DeploymentID) %>%
    arrange(LakeID, deployment_type)
  
  message("\nAvailable sonde deployments:")
  print(head(available_deployments, 20))
  
  message("\nSample-level dataframe created with:")
  message("- ", nrow(matched_data), " zooplankton samples")
  message("- ", length(unique(matched_data$LakeID)), " lakes")
  message("- ", ncol(matched_data), " variables")
  message("- Sonde data: deployment-matched where possible")

  # === Required: Join infection data (by lake + shallow/deep; fallback for single) ===
  if (is.null(infection_path) || !file.exists(infection_path)) {
    stop("Required infection file not found at: ", infection_path)
  } else {
    inf_df <- process_infection_data(infection_path)
    inf_cols <- setdiff(names(inf_df), c("MatchID", "LakeID", "infection_sample_type"))

    # Join exact (LakeID + zoop_sample_type) using a MatchID string key
    joined <- matched_data %>%
      mutate(
        infection_matchid_exact = paste(LakeID, zoop_sample_type, sep = "_"),
        infection_matchid_deep = paste(LakeID, "deep", sep = "_"),
        infection_matchid_shallow = paste(LakeID, "shallow", sep = "_")
      ) %>%
      left_join(inf_df %>% dplyr::select(MatchID, dplyr::all_of(inf_cols)),
                by = c("infection_matchid_exact" = "MatchID"))

    # Add deep + shallow fallbacks (used mainly for zoop_sample_type == 'single')
    inf_deep <- inf_df %>%
      filter(infection_sample_type == "deep") %>%
      dplyr::select(MatchID, dplyr::all_of(inf_cols)) %>%
      rename_with(~ paste0(.x, "_deep"), -MatchID)
    inf_shallow <- inf_df %>%
      filter(infection_sample_type == "shallow") %>%
      dplyr::select(MatchID, dplyr::all_of(inf_cols)) %>%
      rename_with(~ paste0(.x, "_shallow"), -MatchID)

    joined <- joined %>%
      left_join(inf_deep, by = c("infection_matchid_deep" = "MatchID")) %>%
      left_join(inf_shallow, by = c("infection_matchid_shallow" = "MatchID"))

    # Coalesce for SINGLE samples only: exact match (if present) -> deep fallback -> shallow fallback
    is_single <- joined$zoop_sample_type == "single"
    for (col in inf_cols) {
      deep_col <- paste0(col, "_deep")
      shallow_col <- paste0(col, "_shallow")
      if (deep_col %in% names(joined) && shallow_col %in% names(joined)) {
        joined[[col]][is_single] <- dplyr::coalesce(
          joined[[col]][is_single],
          joined[[deep_col]][is_single],
          joined[[shallow_col]][is_single]
        )
      }
    }

    matched_data <- joined %>%
      dplyr::select(
        -infection_matchid_exact, -infection_matchid_deep, -infection_matchid_shallow,
        -dplyr::ends_with("_deep"), -dplyr::ends_with("_shallow")
      )

    message("Infection data joined: added columns: ", paste(inf_cols, collapse = ", "))
  }
  
  # Store additional data as attributes
  attr(matched_data, "sonde_deployments") <- sonde_result
  attr(matched_data, "deployment_summary") <- sonde_result$deployment_summary
  
  return(matched_data)
}