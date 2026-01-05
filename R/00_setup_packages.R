# R/00_setup_packages.R
# --------------------------------------------------------------------
# Load or install CRAN/Bioconductor packages required for:
# - Lipidomics normalization
# - WLCNA (single timepoint)
# - Consensus WLCNA (multi-timepoint)
# - Trait integration + plotting + ML-based downstream functions
# --------------------------------------------------------------------

options(stringsAsFactors = FALSE)

# --------------------------------------------------------------------
# Helper: install if missing, then load
# --------------------------------------------------------------------
install_and_load <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing package: ", pkg)
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# --------------------------------------------------------------------
# CRAN packages
# --------------------------------------------------------------------
cran_packages <- c(
  "caret",
  "randomForest",
  "xgboost",
  "ggplot2",
  "readr",
  "dplyr",
  "writexl",
  "reshape2",
  "pheatmap",
  "viridis",
  "yaml",
  "tibble",
  "stringr"
)

for (pkg in cran_packages) install_and_load(pkg)

# --------------------------------------------------------------------
# Bioconductor packages
# --------------------------------------------------------------------
bioc_packages <- c(
  "WGCNA",
  "ComplexHeatmap",
  "circlize"
)

for (pkg in bioc_packages) install_and_load(pkg, bioc = TRUE)

# --------------------------------------------------------------------
# WGCNA recommended global options
# --------------------------------------------------------------------
allowWGCNAThreads()      # allow multi-threading
enableWGCNAThreads()     # safe activation for most systems

# Preserve numeric precision
options(scipen = 999)

message("All required packages successfully loaded.")
