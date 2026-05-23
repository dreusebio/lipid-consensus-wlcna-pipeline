# R/00b_prepare_raw_data.R
# -------------------------------------------------------
# Prepare raw lipidomics data from the core facility Excel file
#
# Steps:
#   1. Load peak intensity data from the facility Excel
#   2. Load sample-timepoint mapping
#   3. Split 147 samples into 3 timepoints (49 each):
#      Baseline, 36-38wk, Postpartum
#   4. Split lipids into:
#      - Annotated (451 with lipid names) → go into WGCNA pipeline
#      - Unannotated (2,200 without names) → saved separately
#   5. Write raw CSVs in MetaboAnalyst format for each timepoint
#      → ready for 02a_normalize_lipids.R
#
# Input:
#   data/raw/mx_769102_LaSalle_lipidomics_human_dried_blood_spots_12-2023_submit.xlsx
#   data/raw/Seed_Sample_Metabolomics_Subaliquots.xlsx
#
# Output → results/00b_data_preparation/
#   Summary tables, annotation counts, QC checks
# Output → data/raw/raw/ (for 02a_normalize_lipids.R)
#   identified_lipids_baseline_raw.csv
#   identified_lipids_TP36_raw.csv
#   identified_lipids_postpartum_raw.csv
#   unidentified_lipids_baseline_raw.csv
#   unidentified_lipids_TP36_raw.csv
#   unidentified_lipids_postpartum_raw.csv
# -------------------------------------------------------

source("R/00_load_packages.R")
library(openxlsx)

cfg <- yaml::read_yaml("config/config.yml")
out <- cfg$output$s00b
dir.create(out,                       showWarnings = FALSE, recursive = TRUE)
dir.create("data/raw/lipids_raw",     showWarnings = FALSE, recursive = TRUE)
dir.create("data/raw/lipids_normalized", showWarnings = FALSE, recursive = TRUE)

# ── 1. Load raw lipidomics Excel ──────────────────────────
message("Loading lipidomics data from core facility Excel...")
lipid_file <- cfg$data$lipid_excel

wb_data <- read.xlsx(lipid_file, sheet = "Data", colNames = FALSE)

# The Data sheet structure:
#   Rows 1-8:  metadata (Label, Sample#, Species, Organ, etc.)
#   Row 9:     column headers (identifier, annotation, ion species, InChiKey,
#              m/z, ret.time, ESI mode, then sample IDs)
#   Rows 10+:  lipid features

# The Data sheet structure:
#   Row 1 (index 1): "Label", sample_id_1, sample_id_2, ... (e.g. 758160, 758161...)
#   Row 2 (index 2): "Sample #", 1, 2, 3...
#   Rows 3-8:        Species, Organ, Treatment, Sub Treatment, File ID x2
#   Row 9 (index 9): column headers: identifier, annotation, ion species...
#   Rows 10+:        lipid features

# Find the Label row — contains the numeric 6-digit sample IDs
label_row_idx <- NULL
for (i in 1:15) {
  row_vals <- as.character(wb_data[i, ])
  if (any(grepl("^Label$", row_vals, ignore.case = TRUE))) {
    label_row_idx <- i
    break
  }
}
if (is.null(label_row_idx)) stop("Could not find 'Label' row in Data sheet")
message(sprintf("  Found Label row at row index: %d", label_row_idx))

label_row  <- as.character(wb_data[label_row_idx, ])
sample_cols <- which(grepl("^[0-9]{5,7}$", trimws(label_row)))
sample_ids  <- trimws(label_row[sample_cols])

message(sprintf("  Total lipid features: %d", nrow(wb_data) - 9))
message(sprintf("  Total samples found: %d", length(sample_ids)))
message(sprintf("  Sample ID range: %s to %s",
                min(sample_ids), max(sample_ids)))

# Find header row (contains "identifier")
header_row_idx <- NULL
for (i in 1:15) {
  row_vals <- as.character(wb_data[i, ])
  if (any(grepl("^identifier$", row_vals, ignore.case = TRUE))) {
    header_row_idx <- i
    break
  }
}
if (is.null(header_row_idx)) stop("Could not find 'identifier' header row")
message(sprintf("  Found header row at index: %d", header_row_idx))

header <- as.character(wb_data[header_row_idx, ])

# Build lipid data frame from rows after header
lipid_df <- wb_data[(header_row_idx + 1):nrow(wb_data), ]
colnames(lipid_df) <- header
lipid_df <- as.data.frame(lipid_df)

# Rename metadata columns cleanly
colnames(lipid_df)[1:7] <- c("identifier","annotation","ion_species",
                               "InChiKey","mz","ret_time","ESI_mode")

# Rename sample columns to use the 6-digit IDs from the Label row
# The sample columns in lipid_df currently have whatever read.xlsx assigned
# Re-map them using position
lipid_df_sample_cols <- which(colnames(lipid_df) %in% sample_ids |
                               grepl("^[0-9]{5,7}$", colnames(lipid_df)))
if (length(lipid_df_sample_cols) == length(sample_ids)) {
  colnames(lipid_df)[lipid_df_sample_cols] <- sample_ids
} else {
  # Fall back: assign sample_ids by position to columns after col 7
  colnames(lipid_df)[8:(7 + length(sample_ids))] <- sample_ids
}

# Convert intensity columns to numeric
for (sid in sample_ids) {
  if (sid %in% colnames(lipid_df)) {
    lipid_df[[sid]] <- suppressWarnings(as.numeric(lipid_df[[sid]]))
  }
}

message(sprintf("  Lipid data frame: %d features x %d columns",
                nrow(lipid_df), ncol(lipid_df)))

# ── 2. Load timepoint mapping ────────────────────────────
message("\nLoading sample-timepoint mapping...")
tp_map <- read.xlsx(cfg$data$timepoint_map, sheet = 1)
colnames(tp_map) <- c("GlobalID","AliquotID","GROW_ID","Timepoint","AliquotType")
tp_map$AliquotID <- as.character(tp_map$AliquotID)

message(sprintf("  Total samples in map: %d", nrow(tp_map)))
print(table(tp_map$Timepoint))

# Check all sample IDs in lipid file are in the map
missing_from_map <- setdiff(sample_ids, tp_map$AliquotID)
if (length(missing_from_map) > 0) {
  message("  WARNING: ", length(missing_from_map),
          " sample IDs in lipid file not found in timepoint map: ",
          paste(missing_from_map, collapse = ", "))
} else {
  message("  All ", length(sample_ids), " sample IDs matched to timepoint map")
}

# ── 3. Split samples by timepoint ────────────────────────
message("\nSplitting samples by timepoint...")

tp_lookup <- setNames(tp_map$Timepoint, tp_map$AliquotID)
grow_lookup <- setNames(tp_map$GROW_ID, tp_map$AliquotID)

baseline_ids   <- tp_map$AliquotID[tp_map$Timepoint == "Baseline"]
wk36_ids       <- tp_map$AliquotID[tp_map$Timepoint == "36-38wk"]
postpartum_ids <- tp_map$AliquotID[tp_map$Timepoint == "Postpartum"]

# Only keep IDs that exist in the lipid data
baseline_ids   <- intersect(baseline_ids,   sample_ids)
wk36_ids       <- intersect(wk36_ids,       sample_ids)
postpartum_ids <- intersect(postpartum_ids, sample_ids)

message(sprintf("  Baseline:   %d samples", length(baseline_ids)))
message(sprintf("  36-38wk:    %d samples", length(wk36_ids)))
message(sprintf("  Postpartum: %d samples", length(postpartum_ids)))

# ── 4. Split lipids by annotation status ─────────────────
message("\nSplitting lipids by annotation status...")

annotated   <- !is.na(lipid_df$annotation) &
               lipid_df$annotation != "" &
               lipid_df$annotation != " " &
               lipid_df$annotation != "NA"

message(sprintf("  Annotated lipids:   %d", sum(annotated)))
message(sprintf("  Unannotated lipids: %d", sum(!annotated)))

lipid_annotated   <- lipid_df[annotated, ]
lipid_unannotated <- lipid_df[!annotated, ]

# ── 5. Helper: write MetaboAnalyst format CSV ─────────────
# Format:
#   Row 1: "Label", sample_id_1, sample_id_2, ...
#   Row 2: "#CLASS", group_label_1, group_label_2, ...
#   Rows 3+: lipid_name, intensity_1, intensity_2, ...
#
# Group label = GROW participant ID (strips leading zeros for consistency)
write_metaboanalyst_csv <- function(lipid_subset, sample_id_list,
                                    grow_lookup, out_path, tp_label) {
  # Get GROW IDs as group labels
  # Strip "GROW-" prefix and all leading zeros to match trait file rownames
  # e.g. "GROW-00092" → "92", "GROW-00104" → "104"
  group_labels <- grow_lookup[sample_id_list]
  group_labels <- gsub("^GROW-0*", "", group_labels)

  # Feature names = annotation for annotated, identifier for unannotated
  feat_names <- ifelse(
    !is.na(lipid_subset$annotation) & lipid_subset$annotation != "",
    trimws(lipid_subset$annotation),
    lipid_subset$identifier
  )

  # Handle duplicate feature names by appending identifier
  dup <- duplicated(feat_names) | duplicated(feat_names, fromLast = TRUE)
  feat_names[dup] <- paste0(feat_names[dup], "_", lipid_subset$identifier[dup])

  # Extract intensity matrix
  intensity_mat <- lipid_subset[, sample_id_list, drop = FALSE]
  intensity_mat[is.na(intensity_mat)] <- 0  # replace NA with 0 for MetaboAnalyst

  # Build output matrix in MetaboAnalyst format:
  # Row 1: "Label",  sample_id_1, sample_id_2, ...
  # Row 2: "#CLASS", participant_id_1, participant_id_2, ...
  # Rows 3+: lipid_name, intensity_1, intensity_2, ...
  out_mat <- rbind(
    c("Label",  sample_id_list),
    c("#CLASS", group_labels),
    cbind(feat_names, as.matrix(intensity_mat))
  )
  # write.table avoids the extra blank row write.csv adds
  write.table(out_mat, out_path,
              sep = ",", row.names = FALSE, col.names = FALSE,
              quote = TRUE)
  message(sprintf("  Saved %s: %d lipids x %d samples → %s",
                  tp_label, nrow(intensity_mat), length(sample_id_list), out_path))
}

# ── 6. Write annotated (identified) lipid CSVs ───────────
message("\nWriting annotated (identified) lipid CSVs...")

write_metaboanalyst_csv(lipid_annotated, baseline_ids,   grow_lookup,
  "data/raw/lipids_raw/identified_lipids_baseline_raw.csv",   "Identified Baseline")
write_metaboanalyst_csv(lipid_annotated, wk36_ids,       grow_lookup,
  "data/raw/lipids_raw/identified_lipids_TP36_raw.csv",        "Identified 36-38wk")
write_metaboanalyst_csv(lipid_annotated, postpartum_ids, grow_lookup,
  "data/raw/lipids_raw/identified_lipids_postpartum_raw.csv",  "Identified Postpartum")

# ── 7. Write unannotated (unidentified) lipid CSVs ───────
message("\nWriting unannotated (unidentified) lipid CSVs...")

write_metaboanalyst_csv(lipid_unannotated, baseline_ids,   grow_lookup,
  "data/raw/lipids_raw/unidentified_lipids_baseline_raw.csv",   "Unidentified Baseline")
write_metaboanalyst_csv(lipid_unannotated, wk36_ids,       grow_lookup,
  "data/raw/lipids_raw/unidentified_lipids_TP36_raw.csv",        "Unidentified 36-38wk")
write_metaboanalyst_csv(lipid_unannotated, postpartum_ids, grow_lookup,
  "data/raw/lipids_raw/unidentified_lipids_postpartum_raw.csv",  "Unidentified Postpartum")

# ── 8. QC and summary outputs ────────────────────────────
message("\nGenerating QC summary outputs...")

# Annotation summary table
annot_summary <- data.frame(
  Category              = c("Total features", "Annotated", "Unannotated"),
  Count                 = c(nrow(lipid_df), sum(annotated), sum(!annotated)),
  Pct                   = round(c(100,
                                   100*sum(annotated)/nrow(lipid_df),
                                   100*sum(!annotated)/nrow(lipid_df)), 1)
)
write.csv(annot_summary, file.path(out, "annotation_summary.csv"), row.names = FALSE)

# Timepoint summary
tp_summary <- data.frame(
  Timepoint  = c("Baseline","36-38wk","Postpartum"),
  N_samples  = c(length(baseline_ids), length(wk36_ids), length(postpartum_ids))
)
write.csv(tp_summary, file.path(out, "timepoint_summary.csv"), row.names = FALSE)

# Lipid class breakdown for annotated lipids
lipid_annotated$lipid_class <- gsub(" .*","", trimws(lipid_annotated$annotation))
class_counts <- sort(table(lipid_annotated$lipid_class), decreasing = TRUE)
class_df <- data.frame(LipidClass = names(class_counts),
                        Count = as.integer(class_counts))
write.csv(class_df, file.path(out, "annotated_lipid_classes.csv"), row.names = FALSE)

# ESI mode breakdown
esi_counts <- table(lipid_annotated$ESI_mode)
message("\nAnnotated lipids by ESI mode:")
print(esi_counts)
write.csv(as.data.frame(esi_counts), file.path(out, "esi_mode_counts.csv"),
          row.names = FALSE)

# QC plot: lipid class distribution
pdf(file.path(out, "annotated_lipid_class_distribution.pdf"), width = 10, height = 6)
par(mar = c(10, 4, 3, 1))
barplot(head(class_counts, 20),
        main = "Top 20 annotated lipid classes",
        ylab = "Count", las = 2, col = viridis::viridis(20),
        cex.names = 0.8)
dev.off()

# Save full annotated lipid metadata
write.csv(lipid_annotated[, c("identifier","annotation","ion_species",
                               "InChiKey","mz","ret_time","ESI_mode")],
          file.path(out, "annotated_lipids_metadata.csv"), row.names = FALSE)

# Save sample-timepoint mapping used
sample_map_used <- data.frame(
  AliquotID       = c(baseline_ids, wk36_ids, postpartum_ids),
  GROW_ID_raw     = grow_lookup[c(baseline_ids, wk36_ids, postpartum_ids)],
  Participant_ID  = gsub("^GROW-0*", "",
                         grow_lookup[c(baseline_ids, wk36_ids, postpartum_ids)]),
  Timepoint       = c(rep("Baseline",   length(baseline_ids)),
                      rep("36-38wk",    length(wk36_ids)),
                      rep("Postpartum", length(postpartum_ids)))
)
write.csv(sample_map_used, file.path(out, "sample_timepoint_mapping.csv"),
          row.names = FALSE)

message("\n=== Summary ===")
print(annot_summary)
print(tp_summary)
message("\nTop lipid classes:")
print(head(class_df, 10))

# ── Validation: confirm participant IDs match trait file format ──
message("\n=== ID consistency check ===")
all_participant_ids <- unique(gsub("^GROW-0*", "", tp_map$GROW_ID))
message("Participant IDs in lipid data (first 10): ",
        paste(sort(all_participant_ids)[1:10], collapse = ", "))
message("These should match rownames in your trait file.")
message("Trait file example IDs: 92, 104, 114, 117, 120, 122, 129...")
message("If they match, cross-referencing between lipid and trait data will work correctly.")

message("\nScript 00b complete → ", out)
message("Raw CSVs ready in data/raw/lipids_raw/ for 02a_normalize_lipids.R")