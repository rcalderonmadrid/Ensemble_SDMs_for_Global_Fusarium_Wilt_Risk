# ==============================================================================
# Utility Functions for Ensemble Species Distribution Modeling
# ==============================================================================

#' Load Required Packages from Requirements File
#'
#' Reads a requirements.txt file, extracts package names (handling both CRAN
#' and GitHub packages), and loads them into the R session using library().
#'
#' @param file Path to requirements.txt file. Default is "requirements.txt".
#' @return Invisible NULL. Packages are loaded as a side effect.
#' @details
#' The requirements file should contain one package per line:
#' - CRAN packages: just the package name (e.g., "dplyr")
#' - GitHub packages: prefix with "github:" (e.g., "github:user/repo")
#' - Comments (lines starting with #) and empty lines are ignored
#'
#' For GitHub packages, only the repository name (basename) is used when
#' loading with library(), as this is the installed package name.

load_packages <- function(file = "requirements.txt") {
  if (!file.exists(file)) {
    stop("Requirements file '", file, "' not found", call. = FALSE)
  }
  
  # Read and clean package list
  pkgs <- readLines(file, warn = FALSE)
  pkgs <- trimws(pkgs)
  pkgs <- pkgs[nzchar(pkgs) & !grepl("^#", pkgs)]

  # Extract package names (handle github: prefix)
  pkg_names <- vapply(pkgs, function(pkg) {
    if (grepl("^github:", pkg)) {
      basename(gsub("^github:", "", pkg))
    } else {
      pkg
    }
  }, character(1), USE.NAMES = FALSE)

  # Load packages
  invisible(lapply(pkg_names, library, character.only = TRUE))
}

#' Create Directory Structure
#'
#' Creates a base directory and multiple subdirectories for organizing
#' project outputs. If directories already exist, they are not recreated.
#'
#' @param base Base directory path
#' @param dirs Character vector of subdirectory names to create within base
#' @return Invisibly returns a character vector of full paths to all created
#'   (or existing) directories

create_dirs <- function(base, dirs) {
  if (!dir.exists(base)) {
    dir.create(base, recursive = TRUE)
  }
  
  full_paths <- file.path(base, dirs)
  
  for (d in full_paths) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
    }
  }

  invisible(full_paths)
}

#' Environmental Constraint Using One-Class SVM
#'
#' Creates a binary environmental suitability mask using one-class Support
#' Vector Machine (SVM). The SVM learns the environmental envelope from
#' occurrence points and classifies the entire landscape. The output masks
#' environmentally suitable areas (value = NA), while unsuitable areas are
#' highlighted (value = 1) for pseudo-absences sampling.
#'
#' @param occ Data frame with occurrence coordinates. Must have columns "x"
#'   and "y" for longitude and latitude.
#' @param env_layer SpatRaster of environmental predictor variables
#' @param cores Number of CPU cores for parallel processing during spatial
#'   prediction. Default is 1.
#' @return SpatRaster with binary mask where:
#'   - 1 = environmentally unsuitable areas (for pseudo-absence sampling)
#'   - NA = environmentally suitable areas (masked out)
#' @details
#' The function:
#' 1. Extracts environmental values at occurrence locations
#' 2. Trains a one-class SVM with radial kernel (nu = 0.1, allowing 10%
#'    outliers) to learn the environmental envelope
#' 3. Predicts across the full environmental raster (SVM returns: 0 =
#'    unsuitable, 1 = suitable)
#' 4. Recodes predictions for use as a constraint mask:
#'    - SVM class 0 (unsuitable) → 1 (suitable for sampling)
#'    - SVM class 1 (suitable) → NA (exclude from sampling)

env_const <- function(occ, env_layer, cores = 1) {
  # Extract environmental values at occurrence points
  data_coords <- occ[, c("x", "y")]
  env_data <- terra::extract(env_layer, data_coords, ID = FALSE)
  
  # Train one-class SVM (nu = 0.1 allows 10% outliers)
  svm_model <- e1071::svm(
    env_data,
    y = NULL,
    type = "one-classification",
    nu = 0.1,
    scale = TRUE,
    kernel = "radial"
  )
  
  # Predict suitability across landscape (returns 0 = unsuitable, 1 = suitable)
  result <- terra::predict(
    env_layer,
    svm_model,
    na.rm = TRUE,
    cpkgs = "e1071",
    cores = cores,
    wopt = list(memmax = 100, steps = 400)
  )
  
  # Recode for use as sampling constraint: 1 = unsuitable areas, NA = suitable
  result[result == 1] <- NA   # Suitable → NA
  result[result == 0] <- 1    # Unsuitable → 1
  
  return(result)
}

#' Jackknife Test for Variable Importance
#'
#' Performs a leave-one-variable-out analysis to assess the contribution of
#' each predictor variable to model performance. For each variable in the
#' specified set, the function:
#' 1. Removes that variable from the dataset
#' 2. Retrains all specified models without that variable
#' 3. Evaluates model performance (ROC or TSS)
#'
#' Variables whose removal leads to higher performance contribute less unique
#' information and may be redundant.
#'
#' @param data BIOMOD.formated.data object containing species and environmental
#'   data
#' @param models_trained Character vector of model algorithm names to train
#'   (e.g., c("RF", "GBM", "MAXNET"))
#' @param metric Evaluation metric to use: "ROC" or "TSS"
#' @param variables Character vector of variable names to test
#' @param partitions Matrix of cross-validation partitions (logical matrix
#'   with rows = observations, columns = CV folds)
#' @param permut Number of permutations for variable importance calculation in
#'   BIOMOD2. Default is 2.
#' @param nb_cpu Number of CPU cores for parallel model training. Default is 1.
#' @param seed_val Random seed for reproducibility. Default is NULL.
#' @return List with three elements:
#'   - results: Data frame with variable names and mean performance without
#'     each variable
#'   - models_without: List of BIOMOD model objects, one per tested variable
#'   - data_without: List of BIOMOD data objects with each variable removed
#' @details
#' This function is used internally by RemoveCorrVar() to decide which of two
#' correlated variables contributes less and should be removed. The variable
#' whose removal causes the smallest decrease (or largest increase) in
#' performance is considered less important.

JK_test <- function(data,
                    models_trained,
                    metric,
                    variables,
                    partitions,
                    permut = 2,
                    nb_cpu = 1,
                    seed_val = NULL) {
  
  n_vars <- length(variables)
  models_without <- vector("list", length = n_vars)
  data_without <- vector("list", length = n_vars)
  results <- matrix(nrow = n_vars, ncol = 2)
  
  # Set column labels based on metric
  col_labels <- if (metric == "ROC") {
    c("Variable", "Train_ROC_without")
  } else {
    c("Variable", "Train_TSS_without")
  }
  
  # Test each variable by removing it
  for (i in seq_along(variables)) {
    # Create data without current variable
    data_modified <- data
    var_idx <- which(variables[i] == colnames(data@data.env.var))
    data_modified@data.env.var <- data@data.env.var[, -var_idx]
    
    # Configure and train models
    model_opts <- biomod2::bm_ModelingOptions(
      data.type = "binary",
      models = models_trained,
      strategy = "bigboss",
      bm.format = data_modified,
      calib.lines = partitions
    )
    
    jk_model <- biomod2::BIOMOD_Modeling(
      data_modified,
      models = models_trained,
      OPT.user = model_opts,
      CV.strategy = "user.defined",
      CV.user.table = partitions,
      CV.do.full.models = FALSE,
      metric.eval = metric,
      var.import = permut,
      scale.models = FALSE,
      nb.cpu = nb_cpu,
      seed.val = seed_val,
      do.progress = FALSE
    )
    
    # Store results (column 9 contains the evaluation metric)
    results[i, 2] <- mean(biomod2::get_evaluations(jk_model)[, 9])
    models_without[[i]] <- jk_model
    data_without[[i]] <- data_modified
  }
  
  # Format results as data frame
  jk_results <- as.data.frame(results, stringsAsFactors = FALSE)
  colnames(jk_results) <- col_labels
  jk_results[[1]] <- variables
  
  return(list(
    results = jk_results,
    models_without = models_without,
    data_without = data_without
  ))
}

#' Remove Highly Correlated Variables
#'
#' Iteratively identifies and removes highly correlated predictor variables
#' using a data-driven approach that considers variable importance. The
#' algorithm:
#' 1. Computes correlation matrix from background environmental samples
#' 2. Calculates variable importance from trained models (normalized to
#'    percentages)
#' 3. Sorts variables by importance (highest first)
#' 4. For each variable, checks if it's correlated (|r| >= cor_th) with any
#'    other variable
#' 5. When correlations found, uses jackknife test to determine which
#'    correlated variable contributes least
#' 6. Removes the least important variable and updates the correlation matrix
#' 7. Repeats until no correlations exceed the threshold
#'
#' @param model BIOMOD.models.out object from initial model training
#' @param data BIOMOD.formated.data object with species and environmental data
#' @param partitions Matrix of cross-validation partitions
#' @param metric Evaluation metric ("ROC" or "TSS")
#' @param bg_sample Data frame of background environmental samples (no x, y
#'   columns) used to compute correlations
#' @param method Correlation method: "spearman" (default) or "pearson"
#' @param cor_th Correlation threshold for removal. Default is 0.7 (|r| >= 0.7
#'   triggers removal)
#' @param models_trained Character vector of model algorithm names
#' @param permut Number of permutations for variable importance. Default is 2.
#' @param nb_cpu Number of CPU cores for parallel processing. Default is 1.
#' @param seed_val Random seed for reproducibility. Default is NULL.
#' @return List with three elements:
#'   - vars: Character vector of removed variable names
#'   - models_corr_var_removed: Updated BIOMOD model object with final
#'     variable set
#'   - data_corr_var_removed: Updated BIOMOD data object with final variable
#'     set
#' @details
#' Variable importance is normalized within each model run (algorithm, PA
#' replicate, CV fold, permutation) to percentages. The median importance
#' across all runs is used for ranking. This ensures fair comparison across
#' different models and runs.
#'
#' The function prioritizes keeping high-importance variables when choosing
#' which correlated variable to remove.

RemoveCorrVar <- function(model,
                          data,
                          partitions,
                          metric,
                          bg_sample,
                          method = "spearman",
                          cor_th = 0.7,
                          models_trained,
                          permut = 2,
                          nb_cpu = 1,
                          seed_val = NULL) {
  
  # Compute correlation matrix
  cor_matrix <- stats::cor(bg_sample, method = method)
  initial_vars <- model@expl.var.names
  removed_vars <- character(0)
  correlation_removed <- FALSE
  
  # Iterate until no correlations exceed threshold
  while (!correlation_removed) {
    cor_matrix <- as.data.frame(cor_matrix)
    
    # Get variable importance and normalize to percentages
    var_imp <- biomod2::get_variables_importance(model)
    vimp <- data.frame()
    
    # Normalize importance scores within each model run
    for (i in seq_along(unique(var_imp$algo))) {
      for (j in seq_along(unique(var_imp$PA))) {
        for (k in seq_along(unique(var_imp$run))) {
          for (l in seq_along(unique(var_imp$rand))) {
            m <- var_imp %>%
              dplyr::filter(
                .data$algo == unique(var_imp$algo)[i],
                .data$PA == unique(var_imp$PA)[j],
                .data$run == unique(var_imp$run)[k],
                .data$rand == unique(var_imp$rand)[l]
              )
            sum_imp <- sum(m$var.imp)
            m$var.imp <- 100 * m$var.imp / sum_imp
            vimp <- dplyr::bind_rows(vimp, m)
          }
        }
      }
    }
    
    # Summarize importance (median) and sort by descending importance
    scores <- vimp %>%
      dplyr::group_by(.data$expl.var) %>%
      dplyr::summarize(
        Permutation_importance = stats::median(.data$var.imp),
        sd = stats::sd(.data$var.imp),
        .groups = "drop"
      ) %>%
      dplyr::rename(Variable = "expl.var") %>%
      dplyr::arrange(dplyr::desc(.data$Permutation_importance))
    
    vars <- scores$Variable
    discarded_variable <- NULL
    
    # Check each variable for high correlations (starting with most important)
    for (i in seq_along(vars)) {
      coeff <- cor_matrix[vars[i]]
      hcv <- rownames(coeff)[abs(coeff) >= cor_th]
      
      # If correlated variables found, run jackknife to decide which to remove
      if (length(hcv) > 1) {
        jk <- JK_test(
          data,
          models_trained = models_trained,
          metric = metric,
          variables = hcv,
          partitions = partitions,
          permut = permut,
          nb_cpu = nb_cpu,
          seed_val = seed_val
        )
        
        # Remove variable whose exclusion gives best (highest) performance
        # Higher performance without a variable means that variable was
        # contributing less unique information
        best_idx <- which.max(jk$results[, 2])
        model <- jk$models_without[[best_idx]]
        data <- jk$data_without[[best_idx]]
        discarded_variable <- as.character(jk$results$Variable[best_idx])
        
        # Update correlation matrix by removing the discarded variable
        cor_matrix[[discarded_variable]] <- NULL
        cor_matrix <- cor_matrix[
          !(rownames(cor_matrix) == discarded_variable),
        ]
        removed_vars <- c(removed_vars, discarded_variable)
        break
      }
    }
    
    # If no variable removed in this iteration, all remaining correlations
    # are below threshold
    if (is.null(discarded_variable)) {
      correlation_removed <- TRUE
    }
  }
  
  return(list(
    vars = setdiff(removed_vars, initial_vars),
    models_corr_var_removed = model,
    data_corr_var_removed = data
  ))
}

#' Measure Spatial Autocorrelation in Environmental Rasters
#'
#' Estimates the range of spatial autocorrelation in environmental predictor
#' variables using empirical variogram analysis. The range represents the
#' distance at which environmental values become statistically independent,
#' which is critical for designing spatially independent cross-validation
#' blocks.
#'
#' @param env_stack SpatRaster of environmental predictor variables
#' @param num_sample Total number of random points to sample across all layers
#'   for variogram fitting. Default is 5000. These are split into 10 batches.
#' @param seed Random seed for reproducible spatial sampling. Default is 200.
#' @param cores Number of CPU cores for parallel processing. Default is 1.
#' @return List with three elements:
#'   - range_degree: Median autocorrelation range in decimal degrees
#'   - range_km: Median autocorrelation range in kilometers
#'   - range_table: Data frame with layer names, range, and sill values for
#'     each environmental variable
#' @details
#' The function:
#' 1. Randomly samples points from the environmental stack (10 batches for
#'    memory efficiency)
#' 2. Removes duplicate locations
#' 3. Fits an empirical variogram model to each environmental layer using
#'    automap::autofitVariogram()
#' 4. Extracts the range parameter from each variogram (distance to sill)
#' 5. Returns the median range across all layers
#'
#' The range is returned in both kilometers (original unit) and decimal
#' degrees. Conversion assumes WGS84 projection where 1 degree ≈ 111.325 km
#' at the equator. This approximation may be less accurate at high latitudes.
#'
#' The median range is typically used as the minimum block size for spatial
#' cross-validation to ensure training and test sets are spatially independent.

spatial_autocor <- function(env_stack,
                            num_sample = 5000L,
                            seed = 200,
                            cores = 1) {
  
  # Randomly sample points from environmental stack (in 10 batches)
  points <- parallel::mclapply(seq_len(10), function(x) {
    set.seed(seed + x)
    sample <- terra::spatSample(
      env_stack,
      size = num_sample / 10,
      method = "random",
      na.rm = TRUE,
      as.df = TRUE,
      xy = TRUE
    )
    return(sample)
  }, mc.cores = cores)
  
  # Combine batches and remove duplicate locations
  points <- points %>%
    dplyr::bind_rows() %>%
    dplyr::distinct(.data$x, .data$y, .keep_all = TRUE)
  
  points <- terra::vect(
    points,
    geom = c("x", "y"),
    crs = "EPSG:4326"
  )
  
  # Fit variogram models for each environmental layer
  nlayer <- terra::nlyr(env_stack)
  vario_list <- parallel::mclapply(seq_len(nlayer), function(x) {
    pt <- sf::as_Spatial(sf::st_as_sf(points[, x]))
    names(pt) <- "target"
    fit_vario <- automap::autofitVariogram(target ~ 1, pt)
    return(fit_vario)
  }, mc.cores = cores)
  
  # Extract range (3rd column) and sill (2nd column) from row 2 of var_model
  # Row 2 contains the spatial component; row 1 is the nugget
  vario_data <- data.frame(
    layers = names(env_stack),
    range = vapply(vario_list, function(v) v$var_model[2, 3], numeric(1)),
    sill = vapply(vario_list, function(v) v$var_model[2, 2], numeric(1))
  )
  
  vario_data <- vario_data[order(vario_data$range), ]
  the_range <- stats::median(vario_data$range)
  
  # Convert km to degrees (WGS84: 1 degree ≈ 111.325 km at equator)
  the_range_degree <- the_range * 1000 / 111325
  
  return(list(
    range_degree = the_range_degree,
    range_km = the_range,
    range_table = vario_data
  ))
}