# ==============================================================================
# Setup Script: Install Required R Packages
# ==============================================================================
#
# Purpose:
#   This script should be run FIRST before any other analysis scripts.
#   It installs all required R packages specified in requirements.txt.
#
# Usage:
#   Run this script from the project root directory:
#   source("scripts/00_setup.R")
#
#   Or from RStudio: Open this file and click "Source"
#
# Requirements:
#   - Internet connection (for downloading packages)
#   - R version 4.0 or higher recommended
#   - Write permissions to R library directory
#
# What happens:
#   1. Loads the install_requirements() function
#   2. Reads requirements.txt in the project root
#   3. Installs all CRAN and GitHub packages
#   4. Reports success or failure for each package
#
# Troubleshooting:
#   If installation fails:
#   - Check your internet connection
#   - Ensure you have appropriate permissions
#   - Try installing failed packages manually: install.packages("package_name")
#   - For GitHub packages, ensure devtools is working:
#     devtools::session_info()
#
# ==============================================================================

# Load the installation function
source("R/install_requirements.R")

# Install all required packages from requirements.txt
cat("Starting package installation...\n")
failed <- install_requirements()

# Report results
if (length(failed) > 0) {
  cat("\n\u26A0 Some packages failed to install.\n")
  cat("Please review errors above and install manually if needed.\n")
} else {
  cat("\n\u2714 Setup complete! All packages installed successfully.\n")
  cat("You can now run the analysis scripts.\n")
}