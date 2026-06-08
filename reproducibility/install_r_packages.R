# reproducibility/install_r_packages.R
# -------------------------------------------------------
# Run once after pixi install:
#   pixi run install-r-pkgs
# -------------------------------------------------------

options(repos = c(CRAN = "https://cloud.r-project.org"))

lib <- .libPaths()[1]
message("Installing into: ", lib)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = lib)

# Force Bioconductor 3.20 and set all repos explicitly
BiocManager::install(version = "3.20", ask = FALSE)
options(repos = BiocManager::repositories())

message("Installing Bioconductor packages individually...")
bioc_pkgs <- c(
  "Biobase",
  "impute",
  "pcaMethods",
  "siggenes",
  "fgsea",
  "edgeR",
  "BiocParallel",
  "RBGL",
  "WGCNA",
  "ComplexHeatmap",
  "limma",
  "sva",
  "preprocessCore",
  "globaltest",
  "genefilter",
  "AnnotationDbi",
  "annotate",
  "multtest",
  "Rgraphviz",
  "ropls"
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing: ", pkg)
    tryCatch(
      BiocManager::install(pkg, lib = lib, ask = FALSE,
                           update = FALSE, force = TRUE),
      error = function(e) message("FAILED: ", pkg, " - ", e$message)
    )
  } else {
    message("Already installed: ", pkg)
  }
}

# crmn
message("Installing crmn...")
tryCatch(
  install.packages("crmn", lib = lib),
  error = function(e) message("FAILED: crmn - ", e$message)
)

# MetaboAnalystR
message("Installing MetaboAnalystR...")
tryCatch(
  remotes::install_github(
    "xia-lab/MetaboAnalystR",
    build = TRUE,
    build_vignettes = FALSE,
    lib = lib
  ),
  error = function(e) message("FAILED: MetaboAnalystR - ", e$message)
)

message("All done. Now run: pixi run session-info")
