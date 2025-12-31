#' Install Required R Packages from Requirements File
#'
#' This function reads a requirements file and installs all specified
#' packages. It supports both CRAN packages and GitHub repositories.
#' The function is designed to make your project reproducible by
#' ensuring all dependencies are installed before running analysis
#' scripts.
#'
#' @param file Character string. Path to the requirements file.
#'   Default is "requirements.txt".
#'
#' @details
#' The requirements file should list one package per line:
#' - CRAN packages: just the package name (e.g., "dplyr")
#' - GitHub packages: prefix with "github:" followed by repo
#'   (e.g., "github:user/package")
#' - Comments: lines starting with "#" are ignored
#' - Empty lines are ignored
#'
#' The function will:
#' 1. Check if the requirements file exists
#' 2. Install devtools if needed (for GitHub packages)
#' 3. Parse the requirements file
#' 4. Install CRAN packages with dependencies
#' 5. Install GitHub packages with dependencies
#' 6. Report which packages were installed, already present, or
#'    failed
#'
#' @return Invisibly returns a character vector of failed package
#'   names (empty if all succeeded)
#'
#' @examples
#' # Basic usage (assumes requirements.txt in working directory)
#' install_requirements()
#'
#' # Specify a different requirements file
#' install_requirements("R/requirements.txt")
#'
#' # Check if any installations failed
#' failed <- install_requirements()
#' if (length(failed) > 0) {
#'   message(
#'     "Manual intervention needed for: ",
#'     paste(failed, collapse = ", ")
#'   )
#' }

install_requirements <- function(file = "requirements.txt") {

  # Check if requirements file exists
  if (!file.exists(file)) {
    stop("Requirements file '", file, "' not found")
  }

  # Ensure devtools is available for GitHub packages
  if (!requireNamespace("devtools", quietly = TRUE)) {
    message("Installing devtools (required for GitHub packages)...")
    install.packages("devtools", dependencies = TRUE)
  }

  # Read and clean package list
  pkgs <- readLines(file, warn = FALSE)
  pkgs <- trimws(pkgs)
  pkgs <- pkgs[pkgs != ""]
  pkgs <- pkgs[!grepl("^#", pkgs)]
  
  if (length(pkgs) == 0) {
    message("No packages found in requirements file")
    return(invisible(NULL))
  }
  
  # Separate CRAN and GitHub packages
  cran_pkgs <- pkgs[!grepl("^github:", pkgs)]
  github_pkgs <- gsub("^github:", "", pkgs[grepl("^github:", pkgs)])
  
  # Remove devtools from CRAN list if present (already installed)
  cran_pkgs <- cran_pkgs[cran_pkgs != "devtools"]
  
  failed_pkgs <- character()

  # Install CRAN packages
  if (length(cran_pkgs) > 0) {
    message("\n=== Installing CRAN packages ===")
    for (pkg in cran_pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message("Installing CRAN package: ", pkg)
        tryCatch(
          {
            install.packages(pkg, dependencies = TRUE)
          },
          error = function(e) {
            message("✖ Failed to install ", pkg, ": ", e$message)
            failed_pkgs <<- c(failed_pkgs, pkg)
          }
        )
      } else {
        message("✔ ", pkg, " already installed")
      }
    }
  }
  
  # Install GitHub packages
  if (length(github_pkgs) > 0) {
    message("\n=== Installing GitHub packages ===")
    for (repo in github_pkgs) {
      pkg_name <- basename(repo)
      if (!requireNamespace(pkg_name, quietly = TRUE)) {
        message("Installing GitHub package: ", repo)
        tryCatch(
          {
            devtools::install_github(
              repo,
              dependencies = TRUE,
              upgrade = "never"
            )
          },
          error = function(e) {
            message("✖ Failed to install ", repo, ": ", e$message)
            failed_pkgs <<- c(failed_pkgs, repo)
          }
        )
      } else {
        message("✔ ", pkg_name, " already installed")
      }
    }
  }
  
  # Summary
  message("\n=== Installation Summary ===")
  if (length(failed_pkgs) > 0) {
    message("✖ Failed packages: ", paste(failed_pkgs, collapse = ", "))
    warning(
      "Some packages failed to install. ",
      "Please check the errors above."
    )
  } else {
    message("✔ All required packages are installed successfully")
  }

  invisible(failed_pkgs)
}