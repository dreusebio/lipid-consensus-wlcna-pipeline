# R/02_prepare_lipids.R
# -------------------------------------------------------
# Load normalized lipidomics, transpose, build multiExpr_ID
# Outputs → results/02_lipids/
# -------------------------------------------------------

source("R/00_load_packages.R")

cfg <- yaml::read_yaml("config/config.yml")

# Use a 02b-specific output folder for transposed lipid files and logs
out <- if (!is.null(cfg$output$s02b)) {
  cfg$output$s02b
} else {
  file.path(dirname(cfg$output$s02), "02b_lipids")
}

dir.create(out, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$output$processed, showWarnings = FALSE, recursive = TRUE)

demographic_data <- readRDS(file.path(cfg$output$processed, "demographic_data_bmi.rds"))

message("Loading normalized lipidomics data...")

message("Demographic data dimensions: ", paste(dim(demographic_data), collapse = " x "))
message("First demographic rownames: ", paste(head(rownames(demographic_data), 10), collapse = ", "))

# ── Helper: check whether sample IDs are bad ──────────────
bad_sample_ids <- function(ids) {
  ids <- as.character(ids)
  
  is.null(ids) ||
    length(ids) == 0 ||
    any(is.na(ids)) ||
    any(ids == "") ||
    any(ids == "NA") ||
    all(grepl("^NA", ids))
}

# ── Helper: clean sample IDs ──────────────────────────────
clean_ids <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA
  x[x == "NA"] <- NA
  return(x)
}

# ── Read and clean lipidomics file ────────────────────────
# Expected structure:
# Row 1 = Label + aliquot/sample IDs
# Row 2 = #CLASS + participant IDs
# Rows 3+ = lipid names + intensity values

read_and_clean <- function(path, demographic_ids = NULL, timepoint_label = "timepoint") {
  
  message("\nReading lipid file for ", timepoint_label, ": ", path)
  
  raw <- read.csv(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )
  
  message("  Raw dimensions: ", paste(dim(raw), collapse = " x "))
  
  # Candidate IDs
  label_ids       <- clean_ids(unlist(raw[1, -1]))
  participant_ids <- clean_ids(unlist(raw[2, -1]))
  
  lipid_names <- as.character(unlist(raw[3:nrow(raw), 1]))
  lipid_names <- make.unique(trimws(lipid_names))
  
  message("  First row-1 IDs: ", paste(head(label_ids, 10), collapse = ", "))
  message("  First row-2 IDs: ", paste(head(participant_ids, 10), collapse = ", "))
  
  # Prefer row 2 participant IDs if valid
  sample_ids <- participant_ids
  
  # If row 2 IDs are bad, try row 1 IDs
  if (bad_sample_ids(sample_ids)) {
    message("  Row-2 participant IDs look invalid. Trying row-1 label IDs.")
    sample_ids <- label_ids
  }
  
  # If both are bad, but dimensions match phenotype, use demographic rownames
  if (bad_sample_ids(sample_ids)) {
    
    if (!is.null(demographic_ids) && length(demographic_ids) == length(label_ids)) {
      message("  Both row-1 and row-2 IDs look invalid.")
      message("  Using demographic rownames because sample counts match.")
      sample_ids <- as.character(demographic_ids)
    } else {
      stop(
        "Could not determine valid sample IDs for ", timepoint_label,
        ". Check the first two rows of the lipid file."
      )
    }
  }
  
  sample_ids <- make.unique(as.character(sample_ids))
  
  int_mat <- raw[3:nrow(raw), -1, drop = FALSE]
  
  # Convert lipid intensities to numeric
  int_mat <- data.frame(
    lapply(int_mat, function(x) as.numeric(as.character(x))),
    check.names = FALSE
  )
  
  # rows = lipids, columns = participants
  rownames(int_mat) <- lipid_names
  colnames(int_mat) <- sample_ids
  
  # Transpose:
  # rows = participants, columns = lipids
  df_t <- as.data.frame(t(int_mat), check.names = FALSE)
  
  # Clean rownames again after transpose
  rownames(df_t) <- make.unique(trimws(as.character(rownames(df_t))))
  
  message("  Transposed dimensions: ", paste(dim(df_t), collapse = " x "))
  message("  First transposed rownames: ", paste(head(rownames(df_t), 10), collapse = ", "))
  
  # Align to phenotype rownames when possible
  if (!is.null(demographic_ids)) {
    
    demographic_ids <- as.character(demographic_ids)
    overlap <- intersect(rownames(df_t), demographic_ids)
    
    message("  Overlap with demographic IDs: ", length(overlap))
    
    if (length(overlap) == 0) {
      
      if (nrow(df_t) == length(demographic_ids)) {
        message("  No ID overlap, but row counts match.")
        message("  Assigning demographic rownames to lipid data by order.")
        rownames(df_t) <- demographic_ids
      } else {
        stop(
          "No overlap between lipid sample IDs and demographic rownames for ",
          timepoint_label,
          ", and row counts do not match."
        )
      }
      
    } else {
      
      # Keep only samples present in demographic data and order by demographic rownames
      keep_ids <- demographic_ids[demographic_ids %in% rownames(df_t)]
      df_t <- df_t[keep_ids, , drop = FALSE]
      
      message("  After phenotype alignment: ", paste(dim(df_t), collapse = " x "))
      message("  First aligned rownames: ", paste(head(rownames(df_t), 10), collapse = ", "))
    }
  }
  
  return(df_t)
}

# ── Read lipid data at each time point ────────────────────
baseline <- read_and_clean(
  cfg$data$lipids$baseline,
  demographic_ids = rownames(demographic_data),
  timepoint_label = "Baseline"
)

wk36 <- read_and_clean(
  cfg$data$lipids$wk36,
  demographic_ids = rownames(demographic_data),
  timepoint_label = "Wk36_38"
)

postpartum <- read_and_clean(
  cfg$data$lipids$postpartum,
  demographic_ids = rownames(demographic_data),
  timepoint_label = "Postpartum"
)

# ── Final ID checks ───────────────────────────────────────
message("\nFinal sample ID checks:")

message("Baseline first IDs: ", paste(head(rownames(baseline), 10), collapse = ", "))
message("Wk36 first IDs: ", paste(head(rownames(wk36), 10), collapse = ", "))
message("Postpartum first IDs: ", paste(head(rownames(postpartum), 10), collapse = ", "))

message("Baseline overlap with phenotype: ", length(intersect(rownames(baseline), rownames(demographic_data))))
message("Wk36 overlap with phenotype: ", length(intersect(rownames(wk36), rownames(demographic_data))))
message("Postpartum overlap with phenotype: ", length(intersect(rownames(postpartum), rownames(demographic_data))))

if (length(intersect(rownames(baseline), rownames(demographic_data))) == 0) {
  stop("Baseline lipid data has no overlap with demographic rownames.")
}

if (length(intersect(rownames(wk36), rownames(demographic_data))) == 0) {
  stop("Wk36 lipid data has no overlap with demographic rownames.")
}

if (length(intersect(rownames(postpartum), rownames(demographic_data))) == 0) {
  stop("Postpartum lipid data has no overlap with demographic rownames.")
}

# ── Save transposed lipid matrices ────────────────────────
write.csv(
  baseline,
  file.path(out, "lipids_baseline_transposed.csv"),
  row.names = TRUE
)

write.csv(
  wk36,
  file.path(out, "lipids_wk36_transposed.csv"),
  row.names = TRUE
)

write.csv(
  postpartum,
  file.path(out, "lipids_postpartum_transposed.csv"),
  row.names = TRUE
)

# ── Build multiExpr_ID object for WGCNA ───────────────────
multiExpr_ID <- list(
  Identified_Baseline     = list(data = baseline),
  Identified_TP36_38weeks = list(data = wk36),
  Identified_Postpartum   = list(data = postpartum)
)

# Check WGCNA-compatible set structure
exprSize_ID <- checkSets(multiExpr_ID)

# ── Status labels for plotting/annotation ─────────────────
status_labels <- demographic_data %>%
  dplyr::select(egwg, ppwr_e, apo, apo_hdp, apo_gdm, apo_other) %>%
  mutate(
    egwg_status      = ifelse(egwg == 1, "EGWG_YES", "EGWG_NO"),
    ppwr_status      = ifelse(ppwr_e == 1, "PPWR_YES", "PPWR_NO"),
    apo_status       = ifelse(apo == 1, "APO_YES", "APO_NO"),
    apo_hdp_status   = ifelse(apo_hdp == 1, "APO_HDP_YES", "APO_HDP_NO"),
    apo_gdm_status   = ifelse(apo_gdm == 1, "APO_GDM_YES", "APO_GDM_NO"),
    apo_other_status = ifelse(apo_other == 1, "APO_OTHER_YES", "APO_OTHER_NO")
  ) %>%
  dplyr::select(ends_with("_status"))

rownames(status_labels) <- rownames(demographic_data)

# ── Save processed objects ────────────────────────────────
saveRDS(
  multiExpr_ID,
  file.path(cfg$output$processed, "multiExpr_ID.rds")
)

saveRDS(
  exprSize_ID,
  file.path(cfg$output$processed, "exprSize_ID.rds")
)

saveRDS(
  status_labels,
  file.path(cfg$output$processed, "lipid_status_labels.rds")
)

message("Script 02 complete → ", out)
message(
  "  Sets: ", exprSize_ID$nSets,
  " | Samples: ", paste(exprSize_ID$nSamples, collapse = " / "),
  " | Lipids: ", exprSize_ID$nGenes
)