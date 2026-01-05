# 02_prepare_lipid_expression_multiExpr_ID.R
source("R/00_setup_packages.R")

# Load normalized data from MetaboAnalyst (as you do)
identified_lipids_baseline_data_normalized <- readr::read_csv(
  "data_raw/identified_lipids_baseline_data_normalized.csv"
)
identified_lipids_TP36_data_normalized <- readr::read_csv(
  "data_raw/identified_lipids_TP36_data_normalized.csv"
)
identified_postpartum_data_normalized <- readr::read_csv(
  "data_raw/identified_postpartum_data_normalized.csv"
)

# Transpose to samples x lipids and rename columns with first row
transpose_and_clean <- function(df) {
  df_t <- as.data.frame(t(df))
  colnames(df_t) <- df_t[1, ]
  df_t <- df_t[-1, ]
  df_t <- df_t[order(rownames(df_t)), ]
  df_t
}

identified_baseline_transposed  <- transpose_and_clean(identified_lipids_baseline_data_normalized)
identified_TP36_transposed      <- transpose_and_clean(identified_lipids_TP36_data_normalized)
identified_postpartum_transposed <- transpose_and_clean(identified_postpartum_data_normalized)

# Save transposed data (as in your current workflow)
write.csv(identified_baseline_transposed,
          "data_processed/identified_baseline_transposed.csv",
          row.names = TRUE)
write.csv(identified_TP36_transposed,
          "data_processed/identified_TP36_transposed.csv",
          row.names = TRUE)
write.csv(identified_postpartum_transposed,
          "data_processed/identified_postpartum_transposed.csv",
          row.names = TRUE)

# Prepare multi-set expression list for consensus WGCNA:
# (matches your multiExpr_ID idea – one element per timepoint)
multiExpr_ID <- list(
  Identified_Baseline   = list(data = as.data.frame(sapply(identified_baseline_transposed, as.numeric),
                                                    row.names = rownames(identified_baseline_transposed))),
  Identified_TP36_38weeks = list(data = as.data.frame(sapply(identified_TP36_transposed, as.numeric),
                                                      row.names = rownames(identified_TP36_transposed))),
  Identified_Postpartum = list(data = as.data.frame(sapply(identified_postpartum_transposed, as.numeric),
                                                    row.names = rownames(identified_postpartum_transposed)))
)

# Number of samples per set (exprSize_ID in your PDF)
exprSize_ID <- list(
  nSets    = length(multiExpr_ID),
  nSamples = sapply(multiExpr_ID, function(x) nrow(x$data))
)

saveRDS(multiExpr_ID, "data_processed/multiExpr_ID.rds")
saveRDS(exprSize_ID, "data_processed/exprSize_ID.rds")
