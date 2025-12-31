#' Install Required R Packages from Requirements File
#'
#' This function reads a requirements file and installs all specified
#' packages. It supports both CRAN packages and GitHub repositories.
#'
#' @param file Character string. Path to the requirements file.
#'   Default is "requirements.txt".
#'
#' @details
#' The requirements file should list one package per line:
#' - CRAN packages: just the package name (e.g., "dplyr")
#' - GitHub packages: prefix with "github:" followed by repo
#'   (e.g., "github:user/package")
#'
#' @return Invisibly returns a character vector of failed package names
#'   (empty if all succeeded)

install_requirements <- function(file = "requirements.txt") {
  
  # Check if requirements file exists
  if (!file.exists(file)) {
    stop("Requirements file '", file, "' not found", call. = FALSE)
  }
  
  # Ensure devtools is available for GitHub packages
  if (!requireNamespace("devtools", quietly = TRUE)) {
    message("Installing devtools (required for GitHub packages)...")
    install.packages("devtools", dependencies = TRUE)
  }
  
  # Read and clean package list
  pkgs <- readLines(file, warn = FALSE)
  pkgs <- trimws(pkgs)
  pkgs <- pkgs[nzchar(pkgs)]  # More explicit than pkgs != ""
  pkgs <- pkgs[!grepl("^#", pkgs)]
  
  if (length(pkgs) == 0) {
    message("No packages found in requirements file")
    return(invisible(character(0)))
  }
  
  # Separate CRAN and GitHub packages
  is_github <- grepl("^github:", pkgs)
  cran_pkgs <- pkgs[!is_github]
  github_repos <- gsub("^github:", "", pkgs[is_github])
  
  # Remove devtools from CRAN list if present (already installed)
  cran_pkgs <- setdiff(cran_pkgs, "devtools")
  
  failed_pkgs <- character(0)
  
  # Install CRAN packages
  if (length(cran_pkgs) > 0) {
    message("\n=== Installing CRAN packages ===")
    for (pkg in cran_pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message("Installing CRAN package: ", pkg)
        tryCatch(
          {
            install.packages(pkg, dependencies = TRUE)
            message("\u2714 Successfully installed ", pkg)
          },
          error = function(e) {
            message("\u2716 Failed to install ", pkg, ": ", conditionMessage(e))
            failed_pkgs <<- c(failed_pkgs, pkg)
          }
        )
      } else {
        message("\u2714 ", pkg, " already installed")
      }
    }
  }
  
  # Install GitHub packages
  if (length(github_repos) > 0) {
    message("\n=== Installing GitHub packages ===")
    for (repo in github_repos) {
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
            message("\u2714 Successfully installed ", pkg_name)
          },
          error = function(e) {
            message(
              "\u2716 Failed to install ",
              repo,
              ": ",
              conditionMessage(e)
            )
            failed_pkgs <<- c(failed_pkgs, repo)
          }
        )
      } else {
        message("\u2714 ", pkg_name, " already installed")
      }
    }
  }
  
  # Summary
  message("\n=== Installation Summary ===")
  if (length(failed_pkgs) > 0) {
    message("\u2716 Failed packages: ", paste(failed_pkgs, collapse = ", "))
    warning(
      "Some packages failed to install. Please check the errors above.",
      call. = FALSE
    )
  } else {
    message("\u2714 All required packages are installed successfully")
  }
  
  invisible(failed_pkgs)
}