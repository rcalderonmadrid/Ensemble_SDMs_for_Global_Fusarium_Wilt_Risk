# ==============================================================================
# Utility Functions for Ensemble Species Distribution Modeling
# ==============================================================================

# Load Required Packages
# Reads package list from requirements.txt and loads them all
load_packages <- function(file = "requirements.txt") {
  
  if (!file.exists(file)) {
    stop("Requirements file '", file, "' not found")
  }
  
  # Read and clean package list
  pkgs <- readLines(file, warn = FALSE)
  pkgs <- trimws(pkgs)
  pkgs <- pkgs[pkgs != "" & !grepl("^#", pkgs)]
  
  # Extract package names (handle github: prefix)
  pkg_names <- sapply(pkgs, function(pkg) {
    if (grepl("^github:", pkg)) {
      basename(gsub("^github:", "", pkg))
    } else {
      pkg
    }
  }, USE.NAMES = FALSE)
  
  invisible(lapply(pkg_names, library, character.only = TRUE))
}

# Create Directory Structure
# Creates base directory and subdirectories for organizing outputs
# Args: base - base directory path, dirs - vector of subdirectory names
# Returns: invisibly returns full paths to created directories
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

# Environmental Constraint Using One-Class SVM
# Creates binary suitability mask based on environmental envelope
env_const <- function(occ, env_layer, cores) {
  # Extract environmental values at occurrence points
  data <- occ[, c("x", "y")]
  data <- terra::extract(env_layer, data, ID = FALSE)

  # Train one-class SVM (nu=0.1 allows 10% outliers)
  svm_model <- e1071::svm(
    data,
    y = NULL,
    type = "one-classification",
    nu = 0.1,
    scale = TRUE,
    kernel = "radial"
  )

  # Predict and convert to binary mask
  result <- terra::predict(
    env_layer,
    svm_model,
    na.rm = TRUE,
    cpkgs = "e1071",
    cores = cores,
    wopt = list(memmax = 100, steps = 400)
  )
  result[result == 1] <- NA
  result[result == 0] <- 1
  
  return(result)
}

# Jackknife Test for Variable Importance
# Removes each variable one at a time and re-trains models to
# assess impact on performance
JK_test <- function(data,
                    models_trained,
                    metric,
                    variables,
                    partitions,
                    permut = 2,
                    nb_cpu = 1,
                    seed_val = NULL) {
  
  n <- length(variables)
  models_without <- vector("list", length = n)
  data_without <- vector("list", length = n)
  res <- matrix(nrow = n, ncol = 2)

  # Set labels based on metric
  labels <- if (metric == "ROC") {
    c("Variable", "Train_ROC_without")
  } else {
    c("Variable", "Train_TSS_without")
  }

  # Test each variable by removing it
  for (i in seq_along(variables)) {
    data2 <- data
    data2@data.env.var <- data@data.env.var[
      ,
      -which(variables[i] == colnames(data@data.env.var))
    ]

    # Configure and train models without current variable
    model_opts <- bm_ModelingOptions(
      data.type = "binary",
      models = models_trained,
      strategy = "bigboss",
      bm.format = data2,
      calib.lines = partitions
    )

    jk_model <- BIOMOD_Modeling(
      data2,
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
    res[i, 2] <- mean(get_evaluations(jk_model)[, 9])
    models_without[[i]] <- jk_model
    data_without[[i]] <- data2
  }

  # Format results
  jk_test <- as.data.frame(res, stringAsFactor = FALSE)
  colnames(jk_test) <- labels
  jk_test[1] <- variables

  return(list(
    results = jk_test,
    models_without = models_without,
    data_without = data_without
  ))
}

# Remove Highly Correlated Variables
# Iteratively removes correlated variables that contribute least
# to model performance. Uses jackknife test to decide which
# correlated variable to remove at each iteration
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
  removed_vars <- c()
  correlation_removed <- FALSE

  # Iterate until no correlations exceed threshold
  while (correlation_removed == FALSE) {
    cor_matrix <- as.data.frame(cor_matrix)

    # Get variable importance and normalize to percentages
    var_imp <- get_variables_importance(model)
    vimp <- data.frame()

    for (i in seq_along(unique(var_imp$algo))) {
      for (j in seq_along(unique(var_imp$PA))) {
        for (k in seq_along(unique(var_imp$run))) {
          for (l in seq_along(unique(var_imp$rand))) {
            m <- var_imp %>%
              dplyr::filter(
                algo == unique(var_imp$algo)[i],
                PA == unique(var_imp$PA)[j],
                run == unique(var_imp$run)[k],
                rand == unique(var_imp$rand)[l]
              )
            sum_imp <- sum(m$var.imp)
            imp <- 100 * m$var.imp / sum_imp
            m$var.imp <- imp
            vimp <- rbind(vimp, m)
          }
        }
      }
    }
    # Summarize importance (median) and sort by descending importance
    scores <- vimp %>%
      dplyr::group_by(expl.var) %>%
      dplyr::summarize(
        Permutation_importance = median(var.imp),
        sd = sd(var.imp)
      ) %>%
      dplyr::rename(Variable = expl.var) %>%
      dplyr::arrange(dplyr::desc(Permutation_importance))
    vars <- scores$Variable

    discarded_variable <- NULL

    # Check each variable for high correlations
    for (i in seq_along(vars)) {
      coeff <- cor_matrix[vars[i]]
      hcv <- row.names(coeff)[abs(coeff) >= cor_th]

      # If correlated variables found, run jackknife
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
        # Remove variable whose exclusion gives best performance
        index <- which.max(jk$results[, 2])
        model <- jk$models_without[[index]]
        data <- jk$data_without[[index]]
        discarded_variable <- as.character(jk$results$Variable[index])

        # Update correlation matrix
        cor_matrix[discarded_variable] <- NULL
        cor_matrix <- cor_matrix[
          !(row.names(cor_matrix) == discarded_variable),
        ]
        removed_vars <- c(removed_vars, discarded_variable)
        break
      }
    }

    # If no variable removed, all correlations are below threshold
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

# Measure Spatial Autocorrelation in Environmental Rasters
# Estimates spatial autocorrelation range using empirical variograms.
# This determines the distance at which observations become independent
spatial_autocor <- function(env_stack,
                            num_sample = 5000L,
                            seed = 200,
                            cores = 1) {

  # Randomly sample points from environmental stack
  points <- parallel::mclapply(1:10, function(x) {
    set.seed(seed + x)
    sample <- terra::spatSample(
      env_stack,
      num_sample / 10,
      method = "random",
      na.rm = TRUE,
      as.df = TRUE,
      xy = TRUE
    )
    return(sample)
  }, mc.cores = cores)
  # Combine and remove duplicates
  points <- points %>%
    dplyr::bind_rows() %>%
    dplyr::distinct(x, y, .keep_all = TRUE)
  points <- terra::vect(
    points,
    geom = c("x", "y"),
    crs = "+proj=longlat +datum=WGS84"
  )

  # Fit variogram models for each layer
  nlayer <- terra::nlyr(env_stack)
  vario_list <- parallel::mclapply(1:nlayer, function(x) {
    pt <- sf::as_Spatial(sf::st_as_sf(points[, x]))
    names(pt) <- "target"
    fit_vario <- automap::autofitVariogram(target ~ 1, pt)
    return(fit_vario)
  }, mc.cores = cores)
  # Extract range and sill from variogram models
  vario_data <- data.frame(layers = seq_len(nlayer), range = 1, sill = 1)
  for (v in seq_along(vario_list)) {
    vario_data$layers[v] <- names(env_stack)[v]
    vario_data$range[v] <- vario_list[[v]]$var_model[2, 3]
    vario_data$sill[v] <- vario_list[[v]]$var_model[2, 2]
  }

  vario_data <- vario_data[order(vario_data$range), ]
  the_range <- stats::median(vario_data$range)

  # Convert km to degrees (WGS84: 1 degree ≈ 111.325 km at equator)
  the_range_degree <- the_range * 1000 / 111325

  # Return range estimates and full table
  return(list(
    range_degree = the_range_degree,
    range_km = the_range,
    range_table = vario_data
  ))
}