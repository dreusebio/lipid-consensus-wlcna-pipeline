# R/00_load_packages.R
# -------------------------------------------------------
# Load all packages. Installation handled by pixi.
# -------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(readr); library(stringr); library(yaml)
  library(writexl); library(openxlsx)
  library(ggplot2); library(reshape2); library(pheatmap)
  library(viridis); library(patchwork); library(gridExtra)
  library(grid); library(forcats)
  library(WGCNA); library(ComplexHeatmap); library(circlize)
  library(caret); library(randomForest); library(xgboost)
  library(mixOmics); library(pROC)
})
options(stringsAsFactors = FALSE, scipen = 999)
allowWGCNAThreads(); enableWGCNAThreads()
message("All packages loaded.")