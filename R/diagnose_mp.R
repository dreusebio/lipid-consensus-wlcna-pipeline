# R/diagnose_mp.R
# -------------------------------------------------------
# Run this ONCE to inspect the modulePreservation object
# structure so we can fix the column name errors.
# Uses only 10 permutations — takes ~5 minutes locally.
# -------------------------------------------------------

library(WGCNA)
library(yaml)

cfg           <- yaml::read_yaml("config/config.yml")
wgcna_run_tag <- trimws(readLines(
  file.path(cfg$output$processed, "current_wgcna_run_tag.txt")))

multiExpr_ID     <- readRDS(file.path(cfg$output$processed,
                                       "multiExpr_ID.rds"))
consensusMods_ID <- readRDS(file.path(cfg$output$processed,
  paste0("consensusMods_ID_", wgcna_run_tag, ".rds")))

module_colors <- consensusMods_ID$colors

multiData_pair <- list(
  list(data = multiExpr_ID[["Identified_Baseline"]]$data),
  list(data = multiExpr_ID[["Identified_TP36_38weeks"]]$data)
)

mp <- WGCNA::modulePreservation(
  multiData         = multiData_pair,
  multiColor        = list(module_colors, module_colors),
  referenceNetworks = 1,
  nPermutations     = 10,
  randomSeed        = 12345,
  quickCor          = 1,
  verbose           = 0
)

# ── Print full structure ───────────────────────────────────
cat("\n=== Top-level slots in mp ===\n")
print(names(mp))

cat("\n=== mp$preservation slots ===\n")
print(names(mp$preservation))

cat("\n=== mp$preservation$Z structure ===\n")
print(str(mp$preservation$Z, max.level = 3))

cat("\n=== Z[[1]] slot names ===\n")
print(names(mp$preservation$Z[[1]]))

cat("\n=== Z[[1]][[1]] columns ===\n")
print(colnames(mp$preservation$Z[[1]][[1]]))

cat("\n=== Z[[1]][[2]] columns ===\n")
print(colnames(mp$preservation$Z[[1]][[2]]))

cat("\n=== mp$preservation$observed structure ===\n")
print(names(mp$preservation$observed[[1]]))

cat("\n=== observed[[1]][[1]] columns ===\n")
print(colnames(mp$preservation$observed[[1]][[1]]))

cat("\n=== observed[[1]][[2]] columns ===\n")
print(colnames(mp$preservation$observed[[1]][[2]]))

cat("\n=== First few rows of Z[[1]][[2]] ===\n")
print(head(mp$preservation$Z[[1]][[2]]))

cat("\n=== First few rows of observed[[1]][[2]] ===\n")
print(head(mp$preservation$observed[[1]][[2]]))

# Save the object so we can load it in the fixed script
# without rerunning permutations
out <- file.path(cfg$output$s05, wgcna_run_tag)
dir.create(out, showWarnings = FALSE, recursive = TRUE)
saveRDS(mp, file.path(out, "mp_Baseline_vs_Wk36_38_DIAG.rds"))
cat("\nDiagnostic mp object saved to:", out, "\n")
cat("Paste the output above and share it to fix the column names.\n")