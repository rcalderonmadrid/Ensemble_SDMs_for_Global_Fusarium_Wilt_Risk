# ==============================================================================
# Step 1: Remote Sensing Data Modeling - Initial Variable Selection
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
#
# Outputs:
#   - Filtered occurrences with spatial blocks
#   - Pseudo-absence datasets (10 replicates)
#   - Initial BIOMOD models
#   - Variables after correlation removal
#
# ==============================================================================

# Load utility functions and packages
source("R/utils.R")
load_packages()

# ==============================================================================
# Directory Setup
# ==============================================================================
# Create directory structure for data and outputs
data_dirs <- create_dirs(
  base = "Data",
  dirs = c(
    "RS_variables",          # Environmental predictor layers
    "occurrences",           # Species occurrence data
    "step1",                 # Step 1 outputs
    "step1/pseudoabsences"   # Pseudo-absence datasets
  )
)

# Assign directories to variables for easy reference
env_dir <- data_dirs[1]    # Environmental layers
occ_dir <- data_dirs[2]    # Occurrence data
step1_dir <- data_dirs[3]  # Step 1 outputs
psa_dir <- data_dirs[4]    # Pseudo-absences

# ==============================================================================
# 1. Load Environmental Variables
# ==============================================================================
# Load all predictor rasters and stack them into a single SpatRaster object
files_paths <- list.files(
  path = env_dir,
  pattern = "tif$",
  full.names = TRUE
)
env_stack <- terra::rast(files_paths)

# ==============================================================================
# 2. Occurrence Data Preparation
# ==============================================================================
# Read occurrence data (expects columns: x, y for coordinates)
occ <- read.csv(
  paste0(occ_dir, "/occurrences.csv"),
  header = TRUE,
  sep = ";",
  dec = ","
)

# Remove occurrences with missing environmental data
# This ensures all occurrences have complete predictor values
env_dat <- terra::extract(env_stack, occ[, c("x", "y")])
row_has_na <- apply(env_dat, 1, function(x) any(is.na(x)))
occ <- occ[!row_has_na, ]
rm(env_dat, row_has_na)

# ==============================================================================
# 3. Spatial Thinning of Occurrence Records
# ==============================================================================
# Apply geographical filtering to reduce spatial clustering
# This helps address sampling bias and spatial autocorrelation
occ$id <- seq_len(nrow(occ))  # Add unique ID to each occurrence
occ_geofilt <- flexsdm::occfilt_geo(
  data = occ,
  x = "x",
  y = "y",
  env_layer = env_stack,
  method = c("defined", d = "1"),
  prj = terra::crs(env_stack)
)
occ <- occ_geofilt
rm(occ_geofilt)
occ$pr_ab <- 1  # Mark as presences for BIOMOD

# ==============================================================================
# 4. Spatial Autocorrelation Analysis
# ==============================================================================
# Measure spatial autocorrelation to determine appropriate block size for CV.
# Uses variogram analysis to estimate the range of spatial independence
var_autocor <- spatial_autocor(
  env_stack = env_stack,
  num_sample = 500000,
  seed = 200,
  cores = 32
)
saveRDS(var_autocor, paste0(step1_dir, "/spatial_autocorrelation.rds"))

min_block_size_degree <- var_autocor$range_degree
min_block_size_km <- var_autocor$range_km

# ==============================================================================
# 5. Spatial Block Cross-Validation Partitioning
# ==============================================================================
# Create spatial blocks for cross-validation to account for spatial
# autocorrelation. Blocks ensure training/test sets are spatially separated
part_block <- flexsdm::part_sblock(
  env_layer = env_stack,
  data = occ,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  n_part = 5,
  min_res_mult = min_block_size_km,
  max_res_mult = 6000,
  num_grids = 500,
  min_occ = 800,
  prop = 1
)

# Check distribution of occurrences across blocks
part_block$part %>%
  dplyr::group_by(.part) %>%
  dplyr::count()

# Convert block partition to raster matching environmental layers
block_layer <- terra::resample(
  part_block$grid,
  env_stack[[1]],
  method = "near"
)
names(block_layer) <- ".part"
terra::writeRaster(
  block_layer,
  paste0(step1_dir, "/block_partition.tif"),
  overwrite = TRUE
)

# Assign block IDs to occurrences
occ_part <- terra::extract(
  block_layer,
  occ[, c("x", "y")],
  ID = FALSE
)
occ$.part <- occ_part$.part
rm(occ_part)
write.csv(
  occ,
  paste0(step1_dir, "/filtered_occ.csv"),
  row.names = FALSE
)

# ==============================================================================
# 6. Environmental Constraint Layer
# ==============================================================================
# Create environmental constraint using one-class SVM. This masks out
# environmentally unsuitable areas for pseudo-absence sampling
envc_layer <- env_const(occ, env_stack, cores = 25)
terra::writeRaster(
  envc_layer,
  paste0(psa_dir, "/envc_layer.tif"),
  overwrite = TRUE
)

# Create calibration area buffer around occurrences
occ_vect <- terra::vect(
  occ[, c("x", "y")],
  geom = names(occ[, c("x", "y")]),
  crs = terra::crs(env_stack)
)
calibration_area <- terra::buffer(
  occ_vect,
  width = as.numeric(min_block_size_km * 1000)
) %>%
  terra::aggregate()

# Visualize environmental constraint and calibration area
plot(envc_layer)
lines(calibration_area, col = "red", lwd = 5)
points(occ[, c("x", "y")], col = "black", cex = 0.1, pch = 19)

# ==============================================================================
# 7. Pseudo-Absence Generation (10 Replicates)
# ==============================================================================
# Generate multiple pseudo-absence datasets with environmental and geographic
# constraints. Uses parallel processing for efficiency

future::plan(multisession, workers = 2, gc = TRUE)
options(future.globals.maxSize = 25000 * 1024^2)

n_psa_replicates <- 10
psa_rep <- foreach::foreach(
  i = seq_len(n_psa_replicates),
  .options.future = list(seed = TRUE)
) %dofuture% {
  gc()  # Garbage collection to manage memory
  
  # Reload data within parallel worker
  files_paths <- list.files(
    path = env_dir,
    pattern = "tif$",
    full.names = TRUE
  )
  env_stack <- terra::rast(files_paths)
  
  # Create calibration area
  occ_vect <- terra::vect(
    occ[, c("x", "y")],
    geom = names(occ[, c("x", "y")]),
    crs = terra::crs(env_stack)
  )
  ca <- terra::buffer(
    occ_vect,
    width = as.numeric(min_block_size_km * 1000)
  ) %>%
    terra::aggregate()
  rm(occ_vect)
  
  # Crop to calibration area
  block_layer <- terra::rast(paste0(step1_dir, "/block_partition.tif"))
  rlayer <- block_layer %>%
    terra::crop(ca) %>%
    terra::mask(ca)
  rm(ca)
  
  # Geographic constraint: exclude area around occurrences
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
  
  exclusion_radius <- 100 * 1000  # 100 km exclusion radius
  geoc_layer <- geo_const(occ, rlayer, exclusion_radius)
  
  # Combine environmental and geographic constraints
  envc_layer <- terra::rast(paste0(psa_dir, "/envc_layer.tif"))
  
  # Align extents if needed
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
  
  # Create combined constraint layer
  const_layer <- (envc_layer + geoc_layer)
  const_layer <- terra::mask(rlayer, const_layer)
  rm(envc_layer, geoc_layer, rlayer)
  
  # Sample pseudo-absences from constraint layer
  # Sample equal number per spatial block to match presence distribution
  cell_samp <- lapply(seq_along(unique(occ$.part)), function(z) {
    set.seed(150 + i)
    flexsdm::sample_background(
      data = occ,
      x = "x",
      y = "y",
      method = "random",
      n = sum(occ$.part == z),
      rlayer = const_layer,
      maskval = z
    )
  })
  psa_data <- dplyr::bind_rows(cell_samp)
  
  # Assign block IDs to pseudo-absences
  psa_part <- terra::extract(
    block_layer,
    psa_data[, c("x", "y")],
    ID = FALSE
  )
  psa_data$.part <- psa_part$.part
  rm(psa_part, block_layer)
  psa_data$rep <- i
  
  write.csv(
    psa_data,
    paste0(psa_dir, "/psa_", i, ".csv"),
    row.names = FALSE
  )
  rm(const_layer, cell_samp)
  return(psa_data)
}
future::plan(sequential)  # Close parallel session

# Combine all pseudo-absence replicates
psa_rep <- dplyr::bind_rows(psa_rep)

# ==============================================================================
# 8. Prepare Data for BIOMOD2 Modeling
# ==============================================================================
# Format presence and pseudo-absence data according to BIOMOD2 requirements

resp_name <- "Foxy"

# Define dataset dimensions
n_pres <- nrow(occ)
n_psa_per_rep <- nrow(occ)
n_psa_reps <- 10

# Create pseudo-absence table with run assignments
psa_table <- data.frame(
  x = numeric(n_psa_per_rep * n_psa_reps),
  y = numeric(n_psa_per_rep * n_psa_reps),
  stringsAsFactors = FALSE
)

# Add columns for each run
for (run_idx in seq_len(n_psa_reps)) {
  psa_table[[paste0("RUN", run_idx)]] <- FALSE
}

# Fill pseudo-absence table with coordinates for each replicate
for (k in seq_len(n_psa_reps)) {
  row_start <- (k - 1) * n_psa_per_rep + 1
  row_end <- k * n_psa_per_rep
  
  psa_table[row_start:row_end, c("x", "y")] <- 
    psa_rep[row_start:row_end, c("x", "y")]
  psa_table[row_start:row_end, 2 + k] <- TRUE
}

# Create response variable: 1 = presence, NA = pseudo-absence
resp_var <- c(
  rep(1, n_pres),
  rep(NA, n_psa_per_rep * n_psa_reps)
)

# Combine coordinates for presences and pseudo-absences
resp_xy <- rbind(
  occ[, c("x", "y")],
  psa_table[, c("x", "y")]
)

# Create presence table (all runs = TRUE for presences)
pres_table <- data.frame(
  x = occ$x,
  y = occ$y,
  stringsAsFactors = FALSE
)
for (run_idx in seq_len(n_psa_reps)) {
  pres_table[[paste0("RUN", run_idx)]] <- TRUE
}

# Merge presence and pseudo-absence tables
pa_table <- rbind(
  pres_table[, -c(1, 2)],
  psa_table[, -c(1, 2)]
)

# Create BIOMOD2 formatted data object
foxy_biomod_data <- BIOMOD_FormatingData(
  resp.var = resp_var,
  dir.name = step1_dir,
  expl.var = env_stack,
  resp.xy = resp_xy,
  resp.name = resp_name,
  PA.strategy = "user.defined",
  PA.user.table = pa_table,
  na.rm = TRUE
)

saveRDS(
  foxy_biomod_data,
  file = paste0(step1_dir, "/FoxyBiomodData.rds")
)

# ==============================================================================
# 9. Build Cross-Validation Partition Table
# ==============================================================================
# Create partition table matching pseudo-absence replicates.
# Each run uses different pseudo-absences but all presences

n_total_obs <- n_pres + (n_psa_per_rep * n_psa_reps)
partitions <- matrix(NA, nrow = n_total_obs, ncol = n_psa_reps)

for (i in seq_len(n_psa_reps)) {
  # All presences included in every partition
  occ_part <- rep(TRUE, n_pres)
  
  # Only pseudo-absences for this replicate included
  psa_part <- rep(NA, n_psa_per_rep * n_psa_reps)
  psa_indices <- which(psa_rep$rep == i)
  psa_part[psa_indices] <- TRUE
  
  partitions[, i] <- c(occ_part, psa_part)
}

colnames(partitions) <- paste0("_PA", seq_len(n_psa_reps), "_RUN1")
partitions <- as.data.frame(partitions)
partitions[] <- lapply(partitions, as.logical)
partitions <- as.matrix(partitions)

# ==============================================================================
# 10. Initial Model Calibration
# ==============================================================================
# Train multiple algorithm types with default parameters
all_models <- c("ANN", "GBM", "MAXNET", "RF", "GAM", "MARS")

# Configure modeling options (bigboss = default parameters)
model_opts <- bm_ModelingOptions(
  data.type = "binary",
  models = all_models,
  strategy = "bigboss",
  bm.format = foxy_biomod_data,
  calib.lines = partitions
)

# Run model calibration with cross-validation
model_out <- BIOMOD_Modeling(
  foxy_biomod_data,
  modeling.id = "initial",
  models = all_models,
  OPT.user = model_opts,
  CV.strategy = "user.defined",
  CV.user.table = partitions,
  CV.do.full.models = FALSE,
  metric.eval = "TSS",
  var.import = 10,
  scale.models = FALSE,
  nb.cpu = 15,
  seed.val = 150,
  do.progress = TRUE
)

saveRDS(model_out, paste0(step1_dir, "/Model_initial.rds"))

# Summarize model performance
scores <- get_evaluations(model_out) %>%
  dplyr::group_by(algo) %>%
  dplyr::summarise(
    dplyr::across(
      sensitivity:calibration,
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE)
      )
    )
  )

# ==============================================================================
# 11. Correlation-Based Variable Selection
# ==============================================================================
# Generate background sample for computing variable correlations
bg_sample <- parallel::mclapply(
  seq_len(10),
  function(x) {
    set.seed(200 + x)
    sample <- terra::spatSample(
      env_stack,
      100000,
      method = "random",
      na.rm = TRUE,
      as.raster = FALSE,
      as.df = TRUE,
      cells = FALSE,
      xy = TRUE
    )
    return(sample)
  },
  mc.cores = 10
)

# Combine and remove duplicate locations
bg_sample <- bg_sample %>%
  dplyr::bind_rows() %>%
  dplyr::distinct(x, y, .keep_all = TRUE) %>%
  dplyr::select(-c(x, y))

# Remove correlated variables using data-driven approach.
# Iteratively removes variables with |correlation| >= 0.7.
# Uses jackknife test to determine which correlated variable contributes least
rcv <- RemoveCorrVar(
  model = model_out,
  data = foxy_biomod_data,
  partitions = partitions,
  metric = "TSS",
  bg_sample = bg_sample,
  method = "spearman",
  cor_th = 0.7,
  models.trained = all_models,
  permut = 5,
  nb.cpu = 15,
  seed.val = 150
)

# Save final variable set and data
write.csv(
  rcv$models_corr_var_removed@expl.var.names,
  paste0(step1_dir, "/vars_corr_removed.csv"),
  row.names = FALSE
)

saveRDS(
  rcv$data_corr_var_removed,
  paste0(step1_dir, "/data_vars_corr_removed.rds")
)

cat("\n=== Step 1 Complete ===\n")
cat(
  "Final number of variables:",
  length(rcv$models_corr_var_removed@expl.var.names),
  "\n"
)
cat("Removed variables:", paste(rcv$vars, collapse = ", "), "\n")