# ==============================================================================
# Download Optimized Remote Sensing Variables from Figshare
# ==============================================================================
#
# Purpose:
#   Downloads and extracts optimized remote sensing (RS) environmental variables
#   from Figshare repository. These variables have been pre-processed and
#   optimized through Steps 1 and 2 of the SDM workflow.
#
# Requirements:
#   - Internet connection
#   - Sufficient disk space (~500 MB for download + extraction)
#   - Required packages: archive (for .7z extraction)
#
# Outputs:
#   - Data/RS_variables/: Directory containing optimized RS variable rasters
#   - RS_variables.7z: Downloaded archive (can be deleted after)
#
# Note:
#   This script can be used to skip Steps 1-2 if you want to start directly
#   with Step 3 (model tuning and ensemble) using pre-optimized variables.
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Check and Install Required Packages
# ------------------------------------------------------------------------------
cat("Checking required packages...\n")

# Check if archive package is installed (needed for .7z files)
if (!requireNamespace("archive", quietly = TRUE)) {
  cat("Package 'archive' not found. Installing...\n")
  install.packages("archive")
}

# Load archive package
library(archive)

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
# Figshare download URL for optimized RS variables
figshare_url <- "https://figshare.com/ndownloader/files/58967161"

# Define output paths
data_dir <- "Data"
rs_dir <- file.path(data_dir, "RS_variables")
download_file <- file.path(rs_dir, "RS_variables.7z")

# ------------------------------------------------------------------------------
# Create Directory Structure
# ------------------------------------------------------------------------------
cat("\nSetting up directory structure...\n")

# Create RS variables directory if it doesn't exist
if (!dir.exists(rs_dir)) {
  dir.create(rs_dir, recursive = TRUE)
  cat("Created directory:", rs_dir, "\n")
} else {
  cat("Directory already exists:", rs_dir, "\n")
  
  # Check if directory already contains files
  existing_files <- list.files(rs_dir, pattern = "\\.tif$")
  if (length(existing_files) > 0) {
    cat("\nWarning: RS_variables directory contains", length(existing_files), "files.\n")
    cat("Do you want to overwrite? (yes/no): ")
    response <- tolower(trimws(readline()))
    
    if (response != "yes" && response != "y") {
      cat("Download cancelled. Existing files preserved.\n")
      quit(save = "no")
    }
    
    # Remove existing files
    cat("Removing existing files...\n")
    file.remove(file.path(rs_dir, existing_files))
  }
}

# ------------------------------------------------------------------------------
# Download RS Variables from Figshare
# ------------------------------------------------------------------------------
cat("\n=== Downloading Optimized RS Variables ===\n")
cat("Source: Figshare repository\n")
cat("URL:", figshare_url, "\n")
cat("Destination:", download_file, "\n")
cat("File type: 7z compressed archive\n\n")

# Check if file already downloaded
if (file.exists(download_file)) {
  cat("Archive already exists. Skipping download.\n")
} else {
  cat("Starting download... (this may take several minutes)\n")
  
  # Download with progress indication
  tryCatch(
    {
      download.file(
        url = figshare_url,
        destfile = download_file,
        mode = "wb",              # Binary mode for compressed files
        method = "auto",          # Auto-detect best method
        quiet = FALSE             # Show progress
      )
      cat("\nDownload completed successfully!\n")
    },
    error = function(e) {
      cat("\nError downloading file:", conditionMessage(e), "\n")
      cat("Please check your internet connection and try again.\n")
      quit(save = "no", status = 1)
    }
  )
}

# Verify download
if (!file.exists(download_file)) {
  cat("\nError: Downloaded file not found at", download_file, "\n")
  quit(save = "no", status = 1)
}

# Check file size
file_size_mb <- round(file.info(download_file)$size / 1024^2, 2)
cat("Downloaded file size:", file_size_mb, "MB\n")

# ------------------------------------------------------------------------------
# Extract RS Variables from 7z Archive
# ------------------------------------------------------------------------------
cat("\n=== Extracting RS Variables from 7z Archive ===\n")
cat("Extracting to:", rs_dir, "\n")
cat("This may take several minutes...\n\n")

# Extract 7z file using archive package
tryCatch(
  {
    # List contents of archive first (optional, for verification)
    archive_contents <- archive::archive(download_file)
    cat("Archive contains", nrow(archive_contents), "files\n")
    
    # Extract all files
    archive::archive_extract(
      archive = download_file,
      dir = rs_dir
    )
    cat("\nExtraction completed successfully!\n")
  },
  error = function(e) {
    cat("\nError extracting file:", conditionMessage(e), "\n")
    cat("Possible issues:\n")
    cat("  1. The downloaded file may be corrupted\n")
    cat("  2. The 'archive' package may not be properly installed\n")
    cat("  3. Insufficient disk space\n")
    cat("\nTry deleting the .7z file and re-running this script.\n")
    cat("If the problem persists, you may need to extract manually.\n")
    quit(save = "no", status = 1)
  }
)

# ------------------------------------------------------------------------------
# Verify Extraction
# ------------------------------------------------------------------------------
cat("\n=== Verifying Extracted Files ===\n")

# List extracted .tif files
tif_files <- list.files(rs_dir, pattern = "\\.tif$", full.names = FALSE)

if (length(tif_files) == 0) {
  cat("Warning: No .tif files found in", rs_dir, "\n")
  cat("The archive may have an unexpected structure.\n")
  cat("Checking for subdirectories...\n")
  
  # Check for files in subdirectories
  all_files <- list.files(rs_dir, pattern = "\\.tif$", 
                          full.names = FALSE, recursive = TRUE)
  
  if (length(all_files) > 0) {
    cat("Found", length(all_files), ".tif files in subdirectories.\n")
    cat("Moving files to main directory...\n")
    
    # Move files from subdirectories to main directory
    for (f in all_files) {
      file.copy(
        from = file.path(rs_dir, f),
        to = file.path(rs_dir, basename(f)),
        overwrite = TRUE
      )
    }
    
    # Remove subdirectories
    subdirs <- list.dirs(rs_dir, recursive = FALSE)
    unlink(subdirs, recursive = TRUE)
    
    # Update file list
    tif_files <- list.files(rs_dir, pattern = "\\.tif$", full.names = FALSE)
    cat("Successfully moved", length(tif_files), "files to main directory.\n")
  } else {
    cat("No .tif files found. Please check the archive manually.\n")
  }
} else {
  cat("Successfully extracted", length(tif_files), "raster files:\n")
  # Show first 10 files
  display_files <- if (length(tif_files) > 10) {
    c(head(tif_files, 10), paste("... and", length(tif_files) - 10, "more"))
  } else {
    tif_files
  }
  cat(paste0("  - ", display_files, collapse = "\n"), "\n")
}

# Calculate total size of extracted files
if (length(tif_files) > 0) {
  total_size_mb <- sum(
    file.info(file.path(rs_dir, tif_files))$size,
    na.rm = TRUE
  ) / 1024^2
  cat("\nTotal size of extracted files:", round(total_size_mb, 2), "MB\n")
} else {
  total_size_mb <- 0
  cat("\nWarning: No files to calculate size.\n")
}

# ------------------------------------------------------------------------------
# Clean Up
# ------------------------------------------------------------------------------
cat("\n=== Clean Up ===\n")
cat("Do you want to delete the downloaded .7z file to save space? (yes/no): ")
response <- tolower(trimws(readline()))

if (response == "yes" || response == "y") {
  file.remove(download_file)
  cat("Deleted:", download_file, "\n")
  cat("Freed up:", file_size_mb, "MB\n")
} else {
  cat("Keeping download file:", download_file, "\n")
  cat("You can delete it manually later to free up space.\n")
}

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
cat("\n=== Download Summary ===\n")
if (length(tif_files) > 0) {
  cat("✓ RS variables successfully downloaded and extracted\n")
  cat("✓ Location:", rs_dir, "\n")
  cat("✓ Number of raster files:", length(tif_files), "\n")
  cat("✓ Total size:", round(total_size_mb, 2), "MB\n")
  cat("\nYou can now proceed with:\n")
  cat("  - Step 3 (03_RS_step3.R) for model tuning and ensemble\n")
  cat("  - Or Steps 1-2 if you want to customize variable selection\n")
  cat("\n=== Process Complete ===\n")
} else {
  cat("✗ Extraction completed but no .tif files found\n")
  cat("Please check the directory manually:", rs_dir, "\n")
  cat("\n=== Process Completed with Warnings ===\n")
}

