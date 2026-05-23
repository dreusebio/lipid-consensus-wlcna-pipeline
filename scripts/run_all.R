# scripts/run_all.R
# -------------------------------------------------------
# Master pipeline — runs all steps in order
# Usage: pixi run run-all
# Single step: pixi run Rscript R/01_prepare_traits.R
# -------------------------------------------------------
library(yaml)
cfg <- yaml::read_yaml("config/config.yml")
for (d in cfg$output) dir.create(d, showWarnings=FALSE, recursive=TRUE)

steps <- list(
  list(script="R/00b_prepare_raw_data.R",             desc="00b — Split raw data by timepoint & annotation"),
  list(script="R/01_prepare_traits.R",               desc="01 — Prepare trait data"),
  list(script="R/02_normalize_lipids.R",             desc="02a — Normalize raw lipidomics",
       skip_if=!cfg$normalization$run_from_raw),
  list(script="R/02b_prepare_lipids.R",               desc="02 — Prepare lipidomics + build multiExpr"),
  list(script="R/03_run_consensus_wgcna.R",           desc="03 — Run consensus WLCNA"),
  list(script="R/04_module_membership.R",             desc="04 — Module membership & hub lipids"),
  list(script="R/05_module_trait_correlations.R",     desc="05 — Module-trait correlations"),
  list(script="R/06_plot_module_trait_timepoints.R",  desc="06 — Module-trait timepoint plots"),
  list(script="R/07_differential_lipids.R",              desc="07 — QC & differential analysis"),
  list(script="R/09_plsda_permutation.R",             desc="09 — PLS-DA & permutation testing"),
  list(script="R/10_sensitivity_analyses.R",          desc="10 — Sensitivity analyses")
)

start <- Sys.time()
for (i in seq_along(steps)) {
  s    <- steps[[i]]
  skip <- isTRUE(s$skip_if)
  message("\n", strrep("=",60), "\n  Step ", i, "/", length(steps), ": ", s$desc,
          if(skip) " [SKIPPED]", "\n", strrep("=",60))
  if (skip) next
  tryCatch(source(s$script),
           error=function(e) { message("ERROR in ", s$script, ": ", e$message); stop(e) })
}
elapsed <- round(difftime(Sys.time(), start, units="mins"), 1)
message("\nPipeline complete in ", elapsed, " minutes")
sink("reproducibility/session_info.txt"); sessionInfo(); sink()
message("Session info saved.")