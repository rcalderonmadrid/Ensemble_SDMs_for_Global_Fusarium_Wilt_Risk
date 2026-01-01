# ==============================================================================
# Step 1: Correlation-Based Variable Selection
# ==============================================================================
#
# Purpose:
#   Performs initial SDM workflow including:
#   - Data preparation (occurrences, environmental layers)
#   - Spatial filtering and blocking for cross-validation
#   - Pseudo-absence generation with environmental and geographic constraints
#   - Initial model calibration with multiple algorithms
#   - Correlation-based variable selection
#
# Requirements:
#   - Run 00_setup.R first to install required packages
#   - Environmental raster layers in Data/RS_variables/
#     (Note: RS variables not uploaded to Figshare due to large file size.
#      If needed, contact Rocio Calderon at r.calderon@csic.es)
#   - Occurrence data in Data/occurrences/occurrences.csv
#   - Sufficient computational resources (10+ cores recommended)
#
# Outputs:
#   - filtered_occ.csv: Spatially filtered occurrences with block assignments
#   - spatial_autocorrelation.rds: Spatial autocorrelation metrics
#   - block_partition.tif: Spatial blocks for cross-validation
#   - psa_*.csv: Pseudo-absence datasets (10 replicates)
#   - FoxyBiomodData.rds: Formatted biomod2 input object
#   - Model_initial.rds: Initial calibrated models
#   - vars_corr_removed.csv: Variables after correlation filtering
#   - data_vars_corr_removed.rds: Biomod data with filtered variables
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Load Required Libraries and Configure Environment
# ------------------------------------------------------------------------------
# Load utility functions and packages
source("R/utils.R")
load_packages()

# ==============================================================================
# PART 1: DIRECTORY SETUP AND DATA LOADING
# ==============================================================================

# ------------------------------------------------------------------------------
# 1.1 Directory Structure Setup
# ------------------------------------------------------------------------------
# Define paths to existing directories containing input data
env_dir <- "Data/RS_variables"      # Environmental predictor layers
occ_dir <- "Data/occurrences"       # Species occurrence data

# Create directory structure for data and outputs
data_dirs <- create_dirs(
  base = "Data",
  dirs = c(
    "Step 1",                        # Step 1 outputs
    "Step 1/Pseudoabsences"          # Pseudo-absence datasets
  )
)

# Assign directories to variables for easy reference
step1_dir <- data_dirs[1]           # Step 1 outputs
psa_dir <- data_dirs[2]             # Pseudo-absences

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
env_dat <- terra::extract(env_stack, occ[, c("x", "y")])
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
  env_layer = env_stack,
  method = c("defined", d = "1"),   # Defined distance method with 1 km
  prj = terra::crs(env_stack)
)
occ <- occ_geofilt
rm(occ_geofilt)                     # Clean up intermediate object

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
# Uses variogram analysis to estimate the range of spatial independence
var <- spatial_autocor(
  env_stack = env_stack,
  num_sample = 500000,              # Sample 500k points for analysis
  seed = 200,
  cores = 32
)
saveRDS(var, paste0(step1_dir, "/spatial_autocorrelation.rds"))

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
  env_layer = env_stack,
  data = occ,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  n_part = 5,                       # Create 5 spatial blocks (folds)
  min_res_mult = min_block_size_km, # Minimum block size from autocorrelation
  max_res_mult = 6000,              # Maximum block size in km
  num_grids = 500,                  # Number of grid configurations to test
  min_occ = 800,                    # Minimum occurrences per block
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
  env_stack[[1]],                   # Match resolution of first env layer
  method = "near"                   # Nearest neighbor for categorical data
)
names(block_layer) <- ".part"
terra::writeRaster(
  block_layer,
  paste0(step1_dir, "/block_partition.tif"),
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# 3.4 Assign Block IDs to Filtered Occurrences
# ------------------------------------------------------------------------------
# Extract spatial block assignments for each occurrence point
# This links occurrences to their corresponding spatial blocks
occ_blocks <- terra::extract(
  block_layer,
  occ[, c("x", "y")],
  ID = FALSE
)
occ$.part <- occ_blocks$.part
rm(occ_blocks)                      # Clean up temporary variable

# Save filtered occurrences with block assignments
write.csv(
  occ,
  paste0(step1_dir, "/filtered_occ.csv"),
  row.names = FALSE
)

# ==============================================================================
# PART 4: PSEUDO-ABSENCE GENERATION
# ==============================================================================

# ------------------------------------------------------------------------------
# 4.1 Create Environmental Constraint Layer
# ------------------------------------------------------------------------------
# Generate environmental constraint layer to restrict pseudo-absence sampling
# to environmentally unsuitable areas (environmental profiling approach)
# Uses one-class SVM to identify unsuitable environmental conditions
envc_layer <- env_const(occ, env_stack, cores = 25)
terra::writeRaster(
  envc_layer,
  paste0(psa_dir, "/envc_layer.tif"),
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# 4.2 Define Calibration Area
# ------------------------------------------------------------------------------
# Create calibration area buffer around occurrences
# This defines the geographic extent where pseudo-absences can be sampled
occ_vect <- terra::vect(
  occ[, c("x", "y")],
  geom = names(occ[, c("x", "y")]),
  crs = terra::crs(env_stack)
)
calibration_area <- terra::buffer(
  occ_vect,
  width = as.numeric(min_block_size_km * 1000)  # Buffer width in meters
) |>
  terra::aggregate()                # Dissolve overlapping buffers

# Visualize environmental constraint and calibration area
# This helps verify the pseudo-absence sampling region
plot(envc_layer)
lines(calibration_area, col = "red", lwd = 5)
points(occ[, c("x", "y")], col = "black", cex = 0.1, pch = 19)

# ------------------------------------------------------------------------------
# 4.3 Generate Pseudo-Absence Replicates
# ------------------------------------------------------------------------------
# Generate 10 replicate sets of pseudo-absences using parallel processing
# Each replicate uses different random sampling while maintaining:
# - Geographic constraints (minimum distance from presences)
# - Environmental constraints (environmentally unsuitable areas)
# - Spatial block structure (pseudo-absences distributed across blocks)

future::plan(multisession, workers = 2, gc = TRUE)
options(future.globals.maxSize = 25000 * 1024^2)  # 25 GB max object size

psan <- 10                          # Number of pseudo-absence replicates
psa_rep <- foreach::foreach(
  i = seq_len(psan),                # 10 pseudo-absence replicates
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
  
  # Create calibration area buffer around occurrences
  occ_vect <- terra::vect(
    occ[, c("x", "y")],
    geom = names(occ[, c("x", "y")]),
    crs = terra::crs(env_stack)
  )
  ca <- terra::buffer(
    occ_vect,
    width = as.numeric(min_block_size_km * 1000)
  ) |>
    terra::aggregate()
  rm(occ_vect)
  
  # Load and crop block layer to calibration area
  block_layer <- terra::rast(paste0(step1_dir, "/block_partition.tif"))
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
  
  # Load environmental constraint layer
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
    e <- terra::ext(
      apply(df_ext, 1, function(k) k[which.min(abs(k))])
    )
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
  psa <- dplyr::bind_rows(cell_samp)
  
  # Assign block IDs to pseudo-absences
  psa_part <- terra::extract(
    block_layer,
    psa[, c("x", "y")],
    ID = FALSE
  )
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
psa_rep <- dplyr::bind_rows(psa_rep)

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
  x = numeric(psa * psan),          # Initialize x coordinates
  y = numeric(psa * psan),          # Initialize y coordinates
  stringsAsFactors = FALSE
)

# Add columns for each run (replicate)
for (run_idx in seq_len(psan)) {
  psa_table[[paste0("RUN", run_idx)]] <- FALSE
}

# Populate pseudo-absence table with coordinates and replicate assignments
for (k in seq_len(psan)) {
  row_start <- (k - 1) * psa + 1
  row_end <- k * psa
  
  # Assign coordinates for this replicate
  psa_table[row_start:row_end, c("x", "y")] <- 
    psa_rep[row_start:row_end, c("x", "y")]
  # Mark which pseudo-absences belong to this replicate
  psa_table[row_start:row_end, 2 + k] <- TRUE
}

# ------------------------------------------------------------------------------
# 5.3 Create Response Variable Vector
# ------------------------------------------------------------------------------
# Response variable: 1 for presences, NA for pseudo-absences
# biomod2 will populate NAs with 0s based on PA.user.table
resp_var <- c(
  rep(1, pres),
  rep(NA, psa * psan)
)

# ------------------------------------------------------------------------------
# 5.4 Combine Coordinates
# ------------------------------------------------------------------------------
# Combine presence and pseudo-absence coordinates
resp_xy <- rbind(
  occ[, c("x", "y")],
  psa_table[, c("x", "y")]
)

# ------------------------------------------------------------------------------
# 5.5 Create Presence-Absence Table
# ------------------------------------------------------------------------------
# Create presence table (all presences belong to all replicates)
pres_table <- data.frame(
  x = occ$x,
  y = occ$y,
  stringsAsFactors = FALSE
)
for (run_idx in seq_len(psan)) {
  pres_table[[paste0("RUN", run_idx)]] <- TRUE
}

# Merge presence and pseudo-absence tables
# This defines which data points are used in each modeling replicate
pa_table <- rbind(
  pres_table[, -c(1, 2)],
  psa_table[, -c(1, 2)]
)

# ------------------------------------------------------------------------------
# 5.6 Format Data for BIOMOD2
# ------------------------------------------------------------------------------
# Create biomod2 input object with all necessary data
foxy_biomod_data <- BIOMOD_FormatingData(
  resp.var = resp_var,              # Response variable (1/NA)
  dir.name = step1_dir,             # Output directory
  expl.var = env_stack,             # Environmental predictors
  resp.xy = resp_xy,                # Coordinates
  resp.name = resp_name,            # Species name
  PA.strategy = "user.defined",     # User-defined pseudo-absences
  PA.user.table = pa_table,         # PA assignment table
  na.rm = TRUE                      # Remove NAs from predictors
)

# Save biomod2 data object for future use
saveRDS(
  foxy_biomod_data,
  file = paste0(step1_dir, "/FoxyBiomodData.rds")
)

# ==============================================================================
# PART 6: SPATIAL CROSS-VALIDATION SETUP
# ==============================================================================

# ------------------------------------------------------------------------------
# 6.1 Build Cross-Validation Partition Table
# ------------------------------------------------------------------------------
# Create partition table matching pseudo-absence replicates
# Each column represents one pseudo-absence replicate
# TRUE = use for training, NA = not in this replicate
# Note: This is simpler than Step 2 as we don't use spatial blocks for CV here

n_total_obs <- pres + (psa * psan)
partitions <- matrix(NA, nrow = n_total_obs, ncol = psan)

# Build partition table: for each pseudo-absence replicate
for (i in seq_len(psan)) {          # Loop through pseudo-absence replicates
  # All presences included in every partition
  occ_part <- rep(TRUE, pres)
  
  # Only pseudo-absences for this replicate included
  psa_part <- rep(NA, psa * psan)
  psa_indices <- which(psa_rep$rep == i)
  psa_part[psa_indices] <- TRUE
  
  # Combine presence and pseudo-absence partitions
  partitions[, i] <- c(occ_part, psa_part)
}

# Finalize partition table
colnames(partitions) <- paste0("_PA", seq_len(psan), "_RUN1")
partitions <- as.data.frame(partitions)
partitions[] <- lapply(partitions, as.logical)
partitions <- as.matrix(partitions)

# ==============================================================================
# PART 7: MODEL CALIBRATION AND VARIABLE SELECTION
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
model_opts <- bm_ModelingOptions(
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
# Cross-validation uses pseudo-absence replicates
model_out <- BIOMOD_Modeling(
  foxy_biomod_data,
  modeling.id = "initial",          # Model run identifier
  models = all_models,              # Algorithms to use
  OPT.user = model_opts,            # Modeling options
  CV.strategy = "user.defined",     # Use custom cross-validation
  CV.user.table = partitions,       # Pseudo-absence replicate partitions
  CV.do.full.models = FALSE,        # Only run cross-validation models
  metric.eval = "TSS",              # Evaluation metric
  var.import = 10,                  # Number of permutations for var importance
  scale.models = FALSE,             # Don't scale predictions
  nb.cpu = 15,                      # Number of CPU cores to use
  seed.val = 150,                   # Random seed for reproducibility
  do.progress = TRUE                # Show progress bar
)

# Save initial model outputs
saveRDS(model_out, paste0(step1_dir, "/Model_initial.rds"))

# ------------------------------------------------------------------------------
# 7.4 Summarize Model Performance
# ------------------------------------------------------------------------------
# Calculate summary statistics for model evaluation metrics
# Group by algorithm and compute mean and standard deviation
scores <- get_evaluations(model_out) |>
  dplyr::group_by(algo) |>
  dplyr::summarise(
    dplyr::across(
      sensitivity:calibration,
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE)
      )
    ),
    .groups = "drop"
  )

# ==============================================================================
# PART 8: CORRELATION-BASED VARIABLE SELECTION
# ==============================================================================

# ------------------------------------------------------------------------------
# 8.1 Generate Background Sample for Correlation Analysis
# ------------------------------------------------------------------------------
# Generate background sample for computing variable correlations
# Randomly sample approximately 1 million points from environmental stack
# Using multiple iterations with different seeds ensures robust sampling
bg_sample <- parallel::mclapply(
  seq_len(10),
  function(x) {
    set.seed(200 + x)               # Unique seed for each iteration
    sample <- terra::spatSample(
      env_stack,
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
bg_sample <- bg_sample |>
  dplyr::bind_rows() |>
  dplyr::distinct(x, y, .keep_all = TRUE) |>
  dplyr::select(-c(x, y))           # Remove coordinates, keep only env values

# ------------------------------------------------------------------------------
# 8.2 Remove Correlated Variables
# ------------------------------------------------------------------------------
# Remove correlated variables using data-driven approach
# Iteratively removes variables with |correlation| >= 0.7
# Uses jackknife test to determine which correlated variable contributes least
# This reduces multicollinearity while retaining most predictive variables
rcv <- RemoveCorrVar(
  model = model_out,                # Initial models
  data = foxy_biomod_data,          # Biomod data object
  partitions = partitions,          # Cross-validation partitions
  metric = "TSS",                   # Metric to optimize (True Skill Statistic)
  bg_sample = bg_sample,            # Background sample for correlations
  method = "spearman",              # Correlation method (non-parametric)
  cor_th = 0.7,                     # Correlation threshold
  models_trained = all_models,      # Algorithms to evaluate
  permut = 5,                       # Permutations for variable importance
  nb_cpu = 15,                      # Number of CPU cores
  seed_val = 150                    # Random seed for reproducibility
)

# ------------------------------------------------------------------------------
# 8.3 Save Final Variable Set
# ------------------------------------------------------------------------------
# Export final set of variables after correlation filtering
# These will be used in Step 2 for variable optimization
write.csv(
  rcv$models_corr_var_removed@expl.var.names,
  paste0(step1_dir, "/vars_corr_removed.csv"),
  row.names = FALSE
)

# Save biomod data object with filtered variables for potential reuse
saveRDS(
  rcv$data_corr_var_removed,
  paste0(step1_dir, "/data_vars_corr_removed.rds")
)

# ------------------------------------------------------------------------------
# 8.4 Print Summary
# ------------------------------------------------------------------------------
cat("\n=== Step 1 Complete ===\n")
cat(
  "Final number of variables:",
  length(rcv$models_corr_var_removed@expl.var.names),
  "\n"
)
cat("Removed variables:", paste(rcv$vars, collapse = ", "), "\n")

# ==============================================================================
# END OF STEP 1
# ==============================================================================