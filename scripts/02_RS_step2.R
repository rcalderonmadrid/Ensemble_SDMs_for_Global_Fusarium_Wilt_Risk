# ==============================================================================
# Step 2: Remote Sensing Data Modeling - Variable Optimization
# ==============================================================================
#
# Purpose:
#   Performs variable optimization workflow building on Step 1 outputs:
#   - Loads Step 1 uncorrelated remote sensing variables
#   - Applies spatial and environmental filtering to occurrences
#   - Pseudo-absence generation with environmental and geographic constraints
#   - Initial model calibration with multiple algorithms
#   - Optimizes variable set by removing low-contribution predictors
#
# Requirements:
#   - Run 00_setup.R first to install required packages
#   - Run 01_RS_step1.R to obtain uncorrelated variables
#   - Environmental raster layers in RS_variables/ directory
#   - Sufficient computational resources (10+ cores recommended)
#
# Outputs:
#   - filtered_occ.csv: Spatially and environmentally filtered occurrences
#   - spatial_autocorrelation.rds: Spatial autocorrelation metrics
#   - block_partition.tif: Spatial blocks for cross-validation
#   - psa_*.csv: Pseudo-absence datasets (10 replicates)
#   - FoxyBiomodData.rds: Formatted biomod2 input object
#   - Model_initial.rds: Initial calibrated models
#   - vars_optimized.csv: Optimized variable set after filtering
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Load Required Libraries and Configure Environment
# ------------------------------------------------------------------------------
# Load utility functions and packages
source("R/utils.R")
load_packages()

# Configure terra options for large raster processing on HPC
# These settings optimize memory usage and processing for large spatial datasets
# Adjust based on your system's available resources
terra::terraOptions(
  tempdir = "temp",               # Directory for temporary files
  memmax = 100,                   # Maximum memory in GB
  steps = 500,                    # Processing steps for large rasters
  verbose = TRUE                  # Print progress messages
)

# ==============================================================================
# PART 1: DIRECTORY SETUP AND DATA LOADING
# ==============================================================================

# ------------------------------------------------------------------------------
# 1.1 Directory Structure Setup
# ------------------------------------------------------------------------------
# Define paths to existing directories containing input data
env_dir <- "Data/RS_variables"      # Environmental predictor layers
occ_dir <- "Data/occurrences"       # Species occurrence data
step1_dir <- "Data/Step 1"          # Step 1 outputs (uncorrelated variables)

# Create main directory structure for Steps 2
main_dirs <- create_dirs(
  base = "Data",
  dirs = "Step 2"
)

# Assign main directories to variables for easy reference
step2_dir <- main_dirs[1]           # Current step outputs

# Create Step 2 subdirectories for organized output storage
step2_subdirs <- create_dirs(
  base = step2_dir,
  dirs = c(
    "Pseudoabsences"               # Pseudo-absence data storage
  )
)

# Assign subdirectories to variables
psa_dir <- step2_subdirs[1]         # Pseudo-absence directory

# ------------------------------------------------------------------------------
# 1.2 Load Environmental Variables
# ------------------------------------------------------------------------------
# Load all predictor rasters and stack them into a single SpatRaster object
# This creates a multi-layer raster with all environmental variables
files_paths <- list.files(
  path = env_dir,
  pattern = "tif$",                 # Match only .tif files
  full.names = TRUE                 # Return full file paths
)
env_stack <- terra::rast(files_paths)

# Subset to variables that passed correlation filtering in Step 1
# This reduces multicollinearity among predictors
vars_corr_removed <- read.csv(
  paste0(step1_dir, "/vars_corr_removed.csv")
)[, 1]
env_stack_subset <- terra::subset(env_stack, vars_corr_removed)
rm(env_stack)                    # Remove full stack to free memory

# ------------------------------------------------------------------------------
# 1.3 Load and Clean Occurrence Data
# ------------------------------------------------------------------------------
# Read species occurrence data
# Expected format: CSV with 'x' and 'y' columns for longitude and latitude
occ <- read.csv(
  paste0(occ_dir, "/occurrences.csv"),
  header = TRUE,
  sep = ";",
  dec = ","
)

# Remove occurrences with missing environmental data
# This ensures all occurrence points have complete predictor values
# which is essential for model training
env_dat <- terra::extract(env_stack_subset, occ[, c("x", "y")])
row_has_na <- apply(env_dat, 1, function(x) any(is.na(x)))
occ <- occ[!row_has_na, ]
rm(env_dat, row_has_na)             # Clean up temporary variables

# ==============================================================================
# PART 2: OCCURRENCE DATA FILTERING
# ==============================================================================

# ------------------------------------------------------------------------------
# 2.1 Spatial Thinning (Geographic Filtering)
# ------------------------------------------------------------------------------
# Apply geographical filtering to reduce spatial clustering of occurrence points
# This helps address sampling bias and reduces spatial autocorrelation
# Method: maintains minimum distance of 1 km between points
occ$id <- seq_len(nrow(occ))        # Add unique identifier to each occurrence
occ_geofilt <- flexsdm::occfilt_geo(
  data = occ,
  x = "x",
  y = "y",
  env_layer = env_stack_subset,
  method = c("defined", d = "1"),   # Defined distance method with 1 km
  prj = terra::crs(env_stack_subset)
)

# ------------------------------------------------------------------------------
# 2.2 Environmental Filtering
# ------------------------------------------------------------------------------
# Apply environmental filtering to reduce clustering in environmental space
# This ensures occurrences are well-distributed across environmental gradients
# Testing different numbers of bins (2-4) to find optimal filtering

max_bins <- 4
occ_envfilt <- vector("list", length = max_bins - 1)

# Test environmental filtering with 2, 3, and 4 bins
for (i in seq_len(max_bins - 1)) {
  oef <- occfilt_env_edited(
    data = occ_geofilt,
    x = "x",
    y = "y",
    id = "id",
    env_layer = env_stack_subset,
    nbins = i + 1,                  # Number of bins for environmental space
    cores = 10                      # Parallel processing cores
  )
  occ_envfilt[[i]] <- oef
}

# ------------------------------------------------------------------------------
# 2.3 Dimensionality Reduction via PCA
# ------------------------------------------------------------------------------
# Reduce dimensions of environmental space using Principal Component Analysis
# This helps in assessing spatial autocorrelation across environmental gradients

# Randomly sample approximately 1 million points from environmental stack
# Using multiple iterations with different seeds ensures robust sampling
set.seed(150)                       # Set seed for reproducibility
envsample <- parallel::mclapply(
  seq_len(10),
  function(x) {
    set.seed(150 + x)               # Unique seed for each iteration
    sample <- terra::spatSample(
      env_stack_subset,
      size = 100000,                # 100k points per iteration
      method = "random",
      na.rm = TRUE,
      as.raster = FALSE,
      as.df = TRUE,
      cells = FALSE,
      xy = TRUE
    )
    return(sample)
  },
  mc.cores = 10                     # Parallel processing
)

# Combine samples and remove duplicates based on coordinates
envsample <- envsample |>
  dplyr::bind_rows() |>
  dplyr::distinct(x, y, .keep_all = TRUE) |>
  dplyr::select(-c(x, y))           # Remove coordinates, keep only env values

# Perform PCA on environmental sample
# Retain first 3 principal components for spatial autocorrelation analysis
pca <- stats::prcomp(
  envsample,
  retx = TRUE,                      # Return transformed data
  scale. = TRUE,                    # Scale variables to unit variance
  center = TRUE,                    # Center variables to zero mean
  rank. = 3                         # Keep first 3 principal components
)
summary(pca)

# ------------------------------------------------------------------------------
# 2.4 Select Optimal Environmental Filtering Based on Moran's I
# ------------------------------------------------------------------------------
# Select the environmentally filtered occurrences that minimize spatial
# autocorrelation measured by Moran's I index across PC variables
# Lower Moran's I indicates better spatial independence (Velazco et al. 2021)

moran_i <- c()

# Calculate average Moran's I for each filtering scenario
for (i in seq_along(occ_envfilt)) {
  # Extract environmental values at filtered occurrence locations
  pcavar <- terra::extract(
    env_stack_subset,
    occ_envfilt[[i]][c("x", "y")]
  )
  
  # Transform to PC space
  pcaval <- stats::predict(pca, newdata = pcavar)
  
  # Create spatial vector for distance calculations
  coord <- terra::vect(
    occ_envfilt[[i]][c("x", "y")],
    geom = c("x", "y"),
    crs = "+proj=longlat +datum=WGS84"
  )
  
  # Calculate inverse distance matrix for spatial weights
  dist <- terra::distance(coord, coord)
  dist_inv <- 1 / dist
  diag(dist_inv) <- 0               # Set diagonal to zero
  
  # Calculate Moran's I for each of the 3 principal components
  # Average across all PCs to get overall spatial autocorrelation metric
  mi <- mean(
    ape::Moran.I(
      x = pcaval[, 1],
      weight = dist_inv,
      scaled = TRUE
    )$observed,
    ape::Moran.I(
      x = pcaval[, 2],
      weight = dist_inv,
      scaled = TRUE
    )$observed,
    ape::Moran.I(
      x = pcaval[, 3],
      weight = dist_inv,
      scaled = TRUE
    )$observed
  )
  moran_i <- c(moran_i, mi)
}

# Select filtering scenario with minimum Moran's I (best spatial independence)
occ <- occ_envfilt[[which.min(moran_i)]]
rm(occ_geofilt, occ_envfilt)        # Clean up intermediate objects

# Add presence/absence indicator (all filtered points are presences)
occ$pr_ab <- 1

# ==============================================================================
# PART 3: SPATIAL BLOCK PARTITIONING FOR CROSS-VALIDATION
# ==============================================================================

# ------------------------------------------------------------------------------
# 3.1 Measure Spatial Autocorrelation in Environmental Rasters
# ------------------------------------------------------------------------------
# Calculate spatial autocorrelation range in environmental variables
# This determines the minimum block size needed for spatial independence
var <- spatial_autocor(
  env_stack = env_stack_subset,
  num_sample = 500000,              # Sample 500k points for analysis
  seed = 200,
  cores = 25
)
saveRDS(var, paste0(step2_dir, "/spatial_autocorrelation.rds"))

# Extract autocorrelation range for block size determination
min_block_size_degree <- var$range_degree
min_block_size_km <- var$range_km

# ------------------------------------------------------------------------------
# 3.2 Create Spatial Block Partition
# ------------------------------------------------------------------------------
# Generate spatial blocks for cross-validation using flexsdm
# Blocks ensure spatial independence between training and testing data
# This provides more realistic model evaluation than random splits
part_block <- flexsdm::part_sblock(
  env_layer = env_stack_subset,
  data = occ,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  n_part = 5,                       # Create 5 spatial blocks (folds)
  min_res_mult = min_block_size_km, # Minimum block size from autocorrelation
  max_res_mult = 6000,              # Maximum block size in km
  num_grids = 500,                  # Number of grid configurations to test
  min_occ = 810,                    # Minimum occurrences per block
  prop = 1                          # Proportion of data to use
)

# Display distribution of occurrences across blocks
part_block$part |>
  dplyr::group_by(.part) |>
  dplyr::count()

# ------------------------------------------------------------------------------
# 3.3 Convert Block Partition to Raster Layer
# ------------------------------------------------------------------------------
# Transform best block partition to a raster layer matching predictor resolution
# This allows easy extraction of block IDs for any spatial location
block_layer <- terra::resample(
  part_block$grid,
  env_stack_subset[[1]],            # Match resolution of first env layer
  method = "near"                   # Nearest neighbor for categorical data
)
names(block_layer) <- ".part"
terra::writeRaster(
  block_layer,
  paste0(step2_dir, "/block_partition.tif"),
  overwrite = TRUE
)

# ==============================================================================
# PART 4: PSEUDO-ABSENCE GENERATION
# ==============================================================================

# ------------------------------------------------------------------------------
# 4.1 Assign Block IDs to Filtered Occurrences
# ------------------------------------------------------------------------------
# Extract spatial block assignments for each occurrence point
# This links occurrences to their corresponding spatial blocks
occ_blocks <- terra::extract(
  block_layer,
  occ[, c("x", "y")],
  ID = FALSE
)
occ$.part <- occ_blocks$.part
write.csv(
  occ,
  paste0(step2_dir, "/filtered_occ.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 4.2 Define Calibration Area
# ------------------------------------------------------------------------------
# Create calibration area buffer around occurrences
# This defines the geographic extent where pseudo-absences can be sampled
occ_vect <- terra::vect(
  occ[, c("x", "y")],
  geom = c("x", "y"),
  crs = terra::crs(env_stack_subset)
)
calibration_area <- terra::buffer(
  occ_vect,
  width = min_block_size_km * 1000  # Buffer width in meters
) |>
  terra::aggregate()                # Dissolve overlapping buffers

# ------------------------------------------------------------------------------
# 4.3 Create Environmental Constraint Layer
# ------------------------------------------------------------------------------
# Generate environmental constraint layer to restrict pseudo-absence sampling
# to environmentally unsuitable areas (environmental profiling approach)
envc_layer <- env_const(occ, env_stack_subset, cores = 25)
terra::writeRaster(envc_layer, paste0(psa_dir, "/envc_layer.tif"))

# ------------------------------------------------------------------------------
# 4.4 Generate Pseudo-Absence Replicates
# ------------------------------------------------------------------------------
# Generate 10 replicate sets of pseudo-absences using parallel processing
# Each replicate uses different random sampling while maintaining:
# - Geographic constraints (minimum distance from presences)
# - Environmental constraints (environmentally unsuitable areas)
# - Spatial block structure (pseudo-absences distributed across blocks)

future::plan(multisession, workers = 2, gc = TRUE)
options(future.globals.maxSize = 25000 * 1024^2)  # 25 GB max object size

psa_rep <- foreach(
  i = seq_len(10),                  # 10 pseudo-absence replicates
  .options.future = list(seed = TRUE)
) %dofuture% {
  gc()                              # Garbage collection to free memory
  
  # Reload environmental stack within each parallel worker
  files_paths <- list.files(
    path = env_dir,
    pattern = "tif$",
    full.names = TRUE
  )
  env_stack <- terra::rast(files_paths)
  vars_corr_removed <- read.csv(
    paste0(step1_dir, "/vars_corr_removed.csv")
  )[, 1]
  env_stack <- terra::subset(env_stack, vars_corr_removed)
  
  # Define calibration area (4000 km buffer around occurrences)
  occ_vect <- terra::vect(
    occ[, c("x", "y")],
    geom = names(occ[, c("x", "y")]),
    crs = terra::crs(env_stack)
  )
  ca <- terra::buffer(occ_vect, width = as.numeric(4000 * 1000)) |>
    terra::aggregate()
  rm(occ_vect)
  
  # Load and crop block layer to calibration area
  block_layer <- terra::rast(paste0(step2_dir, "/block_partition.tif"))
  rlayer <- block_layer |>
    terra::crop(ca) |>
    terra::mask(ca)
  rm(ca)
  
  # Define geographic constraint function
  # Creates exclusion zones around presence points where pseudo-absences
  # cannot be sampled (minimum distance requirement)
  geo_const <- function(occ, rlayer, exclusion_radius) {
    data <- occ[, c("x", "y")]
    data <- terra::vect(
      data,
      geom = c("x", "y"),
      crs = terra::crs(rlayer)
    )
    b <- terra::buffer(data, width = exclusion_radius)
    b <- terra::rasterize(b, rlayer, background = 0)
    e <- terra::mask(rlayer, b, maskvalues = 1)
    return(e)
  }
  
  # Apply geographic constraint (100 km exclusion radius)
  exclusion_radius <- 100 * 1000    # Convert to meters
  geoc_layer <- geo_const(occ, rlayer, exclusion_radius)
  envc_layer <- terra::rast(paste0(psa_dir, "/envc_layer.tif"))
  
  # Ensure geographic and environmental constraint layers have same extent
  # Crop to common extent if they differ
  if (!all(
    as.vector(terra::ext(geoc_layer)) %in%
      as.vector(terra::ext(envc_layer))
  )) {
    df_ext <- data.frame(
      as.vector(terra::ext(geoc_layer)),
      as.vector(terra::ext(envc_layer))
    )
    e <- terra::ext(apply(df_ext, 1, function(k) k[which.min(abs(k))]))
    geoc_layer <- terra::crop(geoc_layer, e)
    envc_layer <- terra::crop(envc_layer, e)
  }
  
  # Combine environmental and geographic constraints
  # Only cells that pass both constraints can contain pseudo-absences
  const_layer <- (envc_layer + geoc_layer)
  const_layer <- terra::mask(rlayer, const_layer)
  rm(envc_layer, geoc_layer, rlayer)
  
  # Sample pseudo-absences for each spatial block
  # Equal number of pseudo-absences and presences per block
  cell_samp <- lapply(seq_along(unique(occ$.part)), function(z) {
    set.seed(150 + i)               # Ensure reproducibility
    flexsdm::sample_background(
      data = occ,
      x = "x",
      y = "y",
      method = "random",
      n = sum(occ$.part == z),      # Match number of presences in block
      rlayer = const_layer,
      maskval = z                   # Sample only within block z
    )
  })
  
  # Combine pseudo-absences from all blocks
  psa <- cell_samp |> dplyr::bind_rows()
  
  # Assign block IDs to pseudo-absences
  psa_part <- terra::extract(block_layer, psa[c("x", "y")], ID = FALSE)
  psa$.part <- psa_part$.part
  rm(psa_part, block_layer)
  
  # Add replicate identifier
  psa$rep <- i
  
  # Save pseudo-absence replicate
  write.csv(
    psa,
    paste0(psa_dir, "/psa_", i, ".csv"),
    row.names = FALSE
  )
  rm(const_layer, cell_samp)
  return(psa)
}
future::plan(sequential)            # Return to sequential processing

# Combine all pseudo-absence replicates into single data frame
psa_rep <- psa_rep |> dplyr::bind_rows()

# ==============================================================================
# PART 5: BIOMOD2 DATA FORMATTING
# ==============================================================================

# ------------------------------------------------------------------------------
# 5.1 Prepare Response and Predictor Data
# ------------------------------------------------------------------------------
resp_name <- "Foxy"                 # Species/response variable name

# Calculate dimensions
pres <- nrow(occ)                   # Number of presence records
psa <- nrow(occ)                    # Number of pseudo-absences per replicate
psan <- 10                          # Number of pseudo-absence replicates

# ------------------------------------------------------------------------------
# 5.2 Create Pseudo-Absence Table
# ------------------------------------------------------------------------------
# Construct pseudo-absence table for biomod2
# Each column represents one replicate, indicating which pseudo-absences
# belong to that replicate (TRUE/FALSE)
psa_table <- data.frame(
  cbind(
    matrix(1, psa * psan, 2),       # Initialize coordinate columns
    matrix("FALSE", psa * psan, psan)  # Initialize replicate columns
  ),
  stringsAsFactors = FALSE
)
colnames(psa_table) <- c("x", "y", paste("RUN", seq_len(psan), sep = ""))

# Populate pseudo-absence table with coordinates and replicate assignments
start <- 1
for (k in seq_len(psan)) {
  # Assign coordinates for this replicate
  psa_table[seq(start, psa * k), 1:2] <- psa_rep[
    seq(start, psa * k),
    c("x", "y")
  ]
  # Mark which pseudo-absences belong to this replicate
  psa_table[seq(start, psa * k), 2 + k] <- rep("TRUE", psa)
  start <- (psa * k) + 1
}

# Convert coordinates to numeric
psa_table$x <- as.numeric(as.character(psa_table$x))
psa_table$y <- as.numeric(as.character(psa_table$y))

# ------------------------------------------------------------------------------
# 5.3 Create Response Variable Vector
# ------------------------------------------------------------------------------
# Response variable: 1 for presences, NA for pseudo-absences
# biomod2 will populate NAs with 0s based on PA.user.table
resp_var <- as.numeric(c(rep(1, pres), rep(NA, psa * psan)))

# ------------------------------------------------------------------------------
# 5.4 Combine Coordinates
# ------------------------------------------------------------------------------
# Combine presence and pseudo-absence coordinates
resp_xy <- data.frame(rbind(occ[, c("x", "y")], psa_table[, c("x", "y")]))

# ------------------------------------------------------------------------------
# 5.5 Create Presence-Absence Table
# ------------------------------------------------------------------------------
# Create presence table (all presences belong to all replicates)
pres_table <- data.frame(
  cbind(occ[, c("x", "y")]),
  matrix("TRUE", pres, psan)
)
colnames(pres_table) <- c("x", "y", paste("RUN", seq_len(psan), sep = ""))

# Merge presence and pseudo-absence tables
# This defines which data points are used in each modeling replicate
pa_table <- rbind(pres_table[, -c(1, 2)], psa_table[, -c(1, 2)])
pa_table[] <- lapply(pa_table, as.logical)  # Convert to logical

# ------------------------------------------------------------------------------
# 5.6 Format Data for BIOMOD2
# ------------------------------------------------------------------------------
# Create biomod2 input object with all necessary data
foxy_biomod_data <- BIOMOD_FormatingData(
  resp.var = resp_var,              # Response variable
  dir.name = step2_dir,             # Output directory
  expl.var = env_stack_subset,      # Environmental predictors
  resp.xy = resp_xy,                # Coordinates
  resp.name = resp_name,            # Species name
  PA.strategy = "user.defined",     # User-defined pseudo-absences
  PA.user.table = pa_table,         # PA assignment table
  na.rm = TRUE                      # Remove NAs from predictors
)

# Save biomod2 data object for future use
saveRDS(foxy_biomod_data, file = paste0(step2_dir, "/FoxyBiomodData.rds"))

# ==============================================================================
# PART 6: SPATIAL CROSS-VALIDATION SETUP
# ==============================================================================

# ------------------------------------------------------------------------------
# 6.1 Build Cross-Validation Partition Table
# ------------------------------------------------------------------------------
# Create spatial block cross-validation table for biomod2
# Each column represents one fold within one pseudo-absence replicate
# TRUE = use for training, FALSE = use for testing, NA = not in this replicate

parts <- max(unique(occ$.part))     # Number of spatial blocks
partitions <- data.frame(matrix("TRUE", pres + psa * psan, psan * parts))
part_colnames <- c()

# Build partition table: for each replicate and each fold
for (i in seq_len(psan)) {          # Loop through pseudo-absence replicates
  for (j in seq_len(parts)) {       # Loop through spatial folds
    # Presence partition: exclude block j for testing
    occ_part <- ifelse(occ$.part == j, FALSE, TRUE)
    
    # Pseudo-absence partition: exclude block j for testing in replicate i
    psa_part <- matrix("TRUE", psa * psan, 1)
    psa_part[which(psa_rep$rep == i & psa_rep$.part == j), 1] <- FALSE
    psa_part[which(psa_rep$rep != i), 1] <- NA  # NA for other replicates
    
    # Combine presence and pseudo-absence partitions
    part <- c(occ_part, psa_part)
    index <- (i * parts) - (parts - j)
    partitions[, index] <- part
    part_colnames[index] <- paste0("_PA", i, "_RUN", j)
  }
}

# Finalize partition table
colnames(partitions) <- part_colnames
partitions[] <- lapply(partitions, as.logical)
partitions <- as.matrix(partitions)

# ==============================================================================
# PART 7: MODEL CALIBRATION AND VARIABLE OPTIMIZATION
# ==============================================================================

# ------------------------------------------------------------------------------
# 7.1 Define Modeling Algorithms
# ------------------------------------------------------------------------------
# Select algorithms for ensemble modeling
# These represent different modeling approaches and assumptions
all_models <- c(
  "ANN",      # Artificial Neural Networks
  "GBM",      # Generalized Boosted Models
  "MAXNET",   # Maximum Entropy (Maxent)
  "RF",       # Random Forest
  "GAM",      # Generalized Additive Models
  "MARS"      # Multivariate Adaptive Regression Splines
)

# ------------------------------------------------------------------------------
# 7.2 Set Modeling Options
# ------------------------------------------------------------------------------
# Configure modeling parameters using "bigboss" strategy
# This uses recommended default parameters for each algorithm
opt_b <- bm_ModelingOptions(
  data.type = "binary",             # Presence/pseudo-absence data
  models = all_models,
  strategy = "bigboss",             # Use biomod2 default parameters
  bm.format = foxy_biomod_data,
  calib.lines = partitions
)

# ------------------------------------------------------------------------------
# 7.3 Run Initial Model Calibration
# ------------------------------------------------------------------------------
# Calibrate models using all variables and all algorithms
# Cross-validation uses spatial blocks to ensure independence
model_out <- BIOMOD_Modeling(
  foxy_biomod_data,
  modeling.id = "initial",          # Model run identifier
  models = all_models,              # Algorithms to use
  OPT.user = opt_b,                 # Modeling options
  CV.strategy = "user.defined",     # Use custom cross-validation
  CV.user.table = partitions,       # Spatial block partitions
  CV.do.full.models = FALSE,        # Only run cross-validation models
  metric.eval = c("TSS", "ROC"),    # Evaluation metrics
  var.import = 10,                  # Number of permutations for var importance
  scale.models = FALSE,             # Don't scale predictions
  nb.cpu = 15,                      # Number of CPU cores to use
  seed.val = 250,                   # Random seed for reproducibility
  do.progress = TRUE                # Show progress bar
)

# Save initial model outputs
saveRDS(model_out, paste0(step2_dir, "/Model_initial.rds"))

# ------------------------------------------------------------------------------
# 7.4 Variable Optimization
# ------------------------------------------------------------------------------
# Optimize variable set by iteratively removing low-contribution variables
# This process:
# 1. Calculates variable importance for each algorithm
# 2. Removes variables with low importance (below threshold)
# 3. Re-evaluates model performance
# 4. Repeats until all remaining variables are important

varopt <- OptimizeVar(
  model = model_out,                # Initial models
  data = foxy_biomod_data,          # Biomod data object
  partitions = partitions,          # Cross-validation partitions
  metric = "TSS",                   # Metric to optimize (True Skill Statistic)
  th = 2,                           # Importance threshold (%)
  models.trained = c("ANN", "GBM", "MAXNET", "RF", "GAM", "MARS"),
  permut = 10,                      # Permutations for variable importance
  nb.cpu = 15,                      # Number of CPU cores
  seed.val = 100                    # Random seed for reproducibility
)

# ------------------------------------------------------------------------------
# 7.5 Save Optimized Variables
# ------------------------------------------------------------------------------
# Export final set of optimized variables for use in subsequent steps
write.csv(
  varopt$models_var_optimized@expl.var.names,
  paste0(step2_dir, "/vars_optimized.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7.6 Print Summary
# ------------------------------------------------------------------------------
cat("\n=== Step 2 Complete ===\n")
cat(
  "Initial number of variables (from Step 1):",
  length(vars_corr_removed),
  "\n"
)
cat(
  "Final number of optimized variables:",
  length(varopt$models_var_optimized@expl.var.names),
  "\n"
)
cat(
  "Variables removed during optimization:",
  length(vars_corr_removed) - length(varopt$models_var_optimized@expl.var.names),
  "\n"
)
if (length(vars_corr_removed) > length(varopt$models_var_optimized@expl.var.names)) {
  removed_vars <- setdiff(
    vars_corr_removed,
    varopt$models_var_optimized@expl.var.names
  )
  cat("Removed variables:", paste(removed_vars, collapse = ", "), "\n")
}

# ==============================================================================
# END OF STEP 2
# ==============================================================================