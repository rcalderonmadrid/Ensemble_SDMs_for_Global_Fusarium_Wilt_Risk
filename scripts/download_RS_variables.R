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
#   - Extraction tool for .7z:
#       * macOS: p7zip (brew install p7zip)
#       * Linux: p7zip-full / p7zip (via apt/dnf/etc.)
#       * Windows: 7-Zip (can be installed via installr::install.7zip())
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

# Check if httr package is installed (needed for proper HTTP downloads)
if (!requireNamespace("httr", quietly = TRUE)) {
  cat("Package 'httr' not found. Installing...\n")
  install.packages("httr")
}

# Load required packages
library(httr)

# ------------------------------------------------------------------------------
# OS-aware 7z Helper Functions
# ------------------------------------------------------------------------------

get_os <- function() {
  os_type <- tolower(Sys.info()[["sysname"]])
  if (is.na(os_type)) os_type <- .Platform$OS.type
  return(os_type)
}

find_7z <- function() {
  os <- get_os()
  
  # Windows
  if (os == "windows" || .Platform$OS.type == "windows") {
    # Common Windows install locations
    candidates <- c(
      "C:/Program Files/7-Zip/7z.exe",
      "C:/Program Files (x86)/7-Zip/7z.exe",
      file.path(Sys.getenv("ProgramFiles"), "7-Zip", "7z.exe"),
      file.path(Sys.getenv("ProgramFiles(x86)"), "7-Zip", "7z.exe")
    )
    candidates <- candidates[file.exists(candidates)]
    if (length(candidates) > 0) return(candidates[[1]])
    
    # Try PATH
    p <- Sys.which("7z.exe")
    if (nzchar(p)) return(unname(p))
  }
  
  # macOS
  if (os == "darwin") {
    mac_locations <- c(
      "/usr/local/bin/7z",
      "/usr/local/bin/7za",
      "/opt/homebrew/bin/7z",
      "/opt/homebrew/bin/7za"
    )
    for (loc in mac_locations) {
      if (file.exists(loc)) return(loc)
    }
  }
  
  # Linux / Unix (including macOS fallback)
  p <- Sys.which(c("7z", "7za", "7zr"))
  p <- p[nzchar(p)]
  if (length(p) > 0) return(unname(p[[1]]))
  
  return("")
}

install_7z_windows <- function() {
  cat("\n=== Installing 7-Zip on Windows ===\n")
  
  if (!requireNamespace("installr", quietly = TRUE)) {
    cat("Installing 'installr' package...\n")
    install.packages("installr")
  }
  
  cat("Downloading and installing 7-Zip...\n")
  cat("This may take a few minutes and may require administrator privileges...\n\n")
  
  tryCatch({
    installr::install.7zip()
    cat("\n7-Zip installed successfully!\n")
    return(TRUE)
  }, error = function(e) {
    cat("\nAutomatic installation failed:", conditionMessage(e), "\n")
    cat("\nPlease install 7-Zip manually:\n")
    cat("  1. Download from: https://www.7-zip.org/\n")
    cat("  2. Run the installer\n")
    cat("  3. Restart R and run this script again\n")
    return(FALSE)
  })
}

install_7z_mac <- function() {
  cat("\n=== Installing p7zip via Homebrew ===\n")
  
  has_brew <- system("which brew", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  
  if (!has_brew) {
    cat("\nHomebrew is not installed.\n")
    cat("Please install Homebrew first:\n")
    cat('  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"\n\n')
    cat("Then run this script again.\n")
    return(FALSE)
  }
  
  cat("Homebrew detected. Installing p7zip...\n")
  cat("This may take a few minutes...\n\n")
  
  result <- system2("brew", args = c("install", "p7zip"), 
                   stdout = TRUE, stderr = TRUE)
  
  status <- attr(result, "status")
  if (is.null(status) || status == 0) {
    cat("\np7zip installed successfully!\n")
    return(TRUE)
  } else {
    cat("\nInstallation failed. Output:\n")
    cat(paste(result, collapse = "\n"), "\n")
    return(FALSE)
  }
}

install_7z_linux <- function() {
  cat("\n=== Installing p7zip on Linux ===\n")
  
  # Detect package manager
  has_apt <- system("which apt-get", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  has_dnf <- system("which dnf", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  has_yum <- system("which yum", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  has_pac <- system("which pacman", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  has_zyp <- system("which zypper", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  
  cat("Attempting to install p7zip-full...\n")
  cat("This may require sudo password...\n\n")
  
  if (has_apt) {
    cat("Using apt (Debian/Ubuntu)...\n")
    system2("sudo", c("apt-get", "update"), stdout = FALSE, stderr = FALSE)
    result <- system2("sudo", c("apt-get", "install", "-y", "p7zip-full"))
  } else if (has_dnf) {
    cat("Using dnf (Fedora)...\n")
    result <- system2("sudo", c("dnf", "install", "-y", "p7zip", "p7zip-plugins"))
  } else if (has_yum) {
    cat("Using yum (CentOS/RHEL)...\n")
    result <- system2("sudo", c("yum", "install", "-y", "p7zip", "p7zip-plugins"))
  } else if (has_pac) {
    cat("Using pacman (Arch)...\n")
    result <- system2("sudo", c("pacman", "-Sy", "--noconfirm", "p7zip"))
  } else if (has_zyp) {
    cat("Using zypper (openSUSE)...\n")
    result <- system2("sudo", c("zypper", "install", "-y", "p7zip"))
  } else {
    cat("\nNo supported package manager found.\n")
    cat("Please install p7zip manually using your distribution's package manager.\n")
    cat("Common commands:\n")
    cat("  Debian/Ubuntu: sudo apt-get install p7zip-full\n")
    cat("  Fedora: sudo dnf install p7zip p7zip-plugins\n")
    cat("  Arch: sudo pacman -S p7zip\n")
    return(FALSE)
  }
  
  if (result == 0) {
    cat("\np7zip installed successfully!\n")
    return(TRUE)
  } else {
    cat("\nInstallation may have failed (exit code:", result, ")\n")
    return(FALSE)
  }
}

ensure_7z <- function() {
  z <- find_7z()
  if (nzchar(z)) {
    cat("Found 7z at:", z, "\n")
    return(z)
  }
  
  os <- get_os()
  os_name <- switch(os,
                   "windows" = "Windows",
                   "darwin" = "macOS",
                   "linux" = "Linux",
                   "Unknown OS")
  
  cat("\n7z command not found on your system (", os_name, ").\n", sep = "")
  cat("This script requires 7-Zip/p7zip to extract .7z files.\n\n")
  cat("Install 7-Zip/p7zip now? (yes/no): ")
  
  ans <- tolower(trimws(readline()))
  if (!(ans %in% c("y", "yes"))) {
    cat("\nCannot proceed without 7z.\n")
    if (os == "windows") {
      cat("To install manually: https://www.7-zip.org/\n\n")
    } else if (os == "darwin") {
      cat("To install manually: brew install p7zip\n\n")
    } else {
      cat("To install manually: sudo apt-get install p7zip-full (or equivalent)\n\n")
    }
    stop("7z not available")
  }
  
  # Try to install based on OS
  success <- FALSE
  if (os == "windows" || .Platform$OS.type == "windows") {
    success <- install_7z_windows()
  } else if (os == "darwin") {
    success <- install_7z_mac()
  } else if (os == "linux") {
    success <- install_7z_linux()
  } else {
    cat("\nUnsupported operating system for automatic installation.\n")
    cat("Please install 7-Zip/p7zip manually.\n")
  }
  
  if (success) {
    # Wait a moment for installation to complete
    Sys.sleep(2)
    z <- find_7z()
    if (nzchar(z)) {
      cat("\n7z is now available at:", z, "\n")
      return(z)
    }
  }
  
  cat("\nInstallation completed but 7z still not found.\n")
  cat("Please install manually and run this script again.\n")
  stop("7z installation verification failed")
}

extract_with_7z <- function(download_file, rs_dir) {
  cat("\nPreparing to extract archive...\n")
  
  # Ensure 7z is available
  z <- ensure_7z()
  
  cat("Extracting with:", z, "\n")
  cat("From:", download_file, "\n")
  cat("To:", rs_dir, "\n\n")
  
  # Extract with verbose output
  result <- system2(
    z,
    args = c("x", shQuote(download_file), paste0("-o", shQuote(rs_dir)), "-y"),
    stdout = TRUE,
    stderr = TRUE
  )
  
  status <- attr(result, "status")
  
  if (!is.null(status) && status != 0) {
    cat("\nExtraction failed. Output:\n")
    cat(paste(result, collapse = "\n"), "\n")
    stop("7z extraction failed with status: ", status)
  }
  
  cat("\n7z extraction completed!\n")
  invisible(TRUE)
}

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
  cat("Attempting direct download from Figshare...\n\n")
  
  # Download using httr with proper headers - direct attempt
  download_result <- tryCatch(
    {
      response <- httr::GET(
        figshare_url,
        httr::add_headers(
          "User-Agent" = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36",
          "Accept" = "*/*",
          "Accept-Language" = "en-US,en;q=0.9",
          "Referer" = "https://figshare.com/"
        ),
        httr::write_disk(download_file, overwrite = TRUE),
        httr::progress(),
        httr::timeout(3600)  # 1 hour timeout for large files
      )
      
      status <- httr::status_code(response)
      
      if (status == 200) {
        "success"
      } else if (status == 403) {
        list(error = TRUE, 
             message = paste("Access forbidden (403). Figshare may be blocking automated downloads."),
             status = status)
      } else if (status == 404) {
        list(error = TRUE, 
             message = paste("File not found (404). The URL may be incorrect or file was removed."),
             status = status)
      } else {
        list(error = TRUE, 
             message = paste("HTTP status:", status),
             status = status)
      }
    },
    error = function(e) {
      return(list(error = TRUE, message = conditionMessage(e)))
    }
  )
  
  # Check download result
  if (is.list(download_result) && !is.null(download_result$error)) {
    cat("\nError downloading file:", download_result$message, "\n")
    
    if (!is.null(download_result$status) && download_result$status == 403) {
      cat("\nFigshare is blocking automated downloads.\n")
      cat("\nPlease download manually:\n")
      cat("  1. Open this URL in your web browser:\n")
      cat("     ", figshare_url, "\n")
      cat("  2. The download should start automatically\n")
      cat("  3. Save the file as 'RS_variables.7z'\n")
      cat("  4. Move it to:", rs_dir, "\n")
      cat("  5. Run this script again - it will skip download and extract the file\n")
    } else {
      cat("\nAlternative download methods:\n")
      cat("  Option 1 - Manual download:\n")
      cat("     Open in browser:", figshare_url, "\n")
      cat("     Save as 'RS_variables.7z' in:", rs_dir, "\n")
      cat("\n  Option 2 - Command line (if you have curl/wget):\n")
      cat("     cd", rs_dir, "\n")
      cat("     curl -L -o RS_variables.7z", figshare_url, "\n")
      cat("\n  Then run this script again to extract.\n")
    }
    
    # Remove partial download if exists
    if (file.exists(download_file) && file.info(download_file)$size < 1000000) {
      file.remove(download_file)
      cat("\nRemoved incomplete download.\n")
    }
    
    stop("Download failed - please download manually as described above")
  } else {
    cat("\nDownload completed successfully!\n")
  }
}

# Verify download
if (!file.exists(download_file)) {
  cat("\nError: Downloaded file not found at", download_file, "\n")
  cat("Download may have failed or file was saved to a different location.\n")
  stop("Download verification failed")
}

# Check file size
file_size_mb <- round(file.info(download_file)$size / 1024^2, 2)
cat("Downloaded file size:", file_size_mb, "MB\n")

# Verify file is not empty or too small
if (file_size_mb < 1) {
  cat("\nWarning: Downloaded file is suspiciously small (", file_size_mb, "MB)\n")
  cat("The file may be corrupted or the download was incomplete.\n")
  cat("Please try downloading manually from:", figshare_url, "\n")
  stop("Downloaded file appears to be invalid")
}

# Check file signature (7z files start with "7z¼¯'" or hex 37 7A BC AF 27 1C)
cat("\nVerifying file format...\n")
con <- file(download_file, "rb")
file_header <- readBin(con, "raw", n = 6)
close(con)

is_7z <- identical(file_header[1:2], as.raw(c(0x37, 0x7A)))

if (!is_7z) {
  cat("Warning: File does not appear to be a valid 7z archive.\n")
  cat("File header:", paste(as.character(file_header[1:6]), collapse = " "), "\n")
  cat("\nThe downloaded file may be:\n")
  cat("  - An HTML error page instead of the actual file\n")
  cat("  - Corrupted during download\n")
  cat("  - In a different format than expected\n")
  cat("\nPlease download manually from:", figshare_url, "\n")
  cat("Save as 'RS_variables.7z' in:", rs_dir, "\n")
  stop("Invalid file format")
}

cat("File format verified: Valid 7z archive\n")

# ------------------------------------------------------------------------------
# Extract RS Variables from 7z Archive
# ------------------------------------------------------------------------------
cat("\n=== Extracting RS Variables from 7z Archive ===\n")
cat("Extracting to:", rs_dir, "\n")
cat("This may take several minutes...\n\n")

# Always use system 7z/7za (with OS-aware install prompt if missing)
extract_with_7z(download_file, rs_dir)
cat("\nExtraction completed successfully (7z)!\n")

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

# Safely read user input
response <- tryCatch({
  tolower(trimws(readline()))
}, error = function(e) {
  "no"  # Default to keeping file if there's an error
})

# Handle empty input
if (is.null(response) || response == "") {
  response <- "no"
}

if (response %in% c("yes", "y")) {
  tryCatch({
    file.remove(download_file)
    cat("Deleted:", download_file, "\n")
    cat("Freed up:", file_size_mb, "MB\n")
  }, error = function(e) {
    cat("Could not delete file:", conditionMessage(e), "\n")
    cat("You can delete it manually:", download_file, "\n")
  })
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