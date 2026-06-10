# R/03_normalize_lipids.R
# -------------------------------------------------------
# Normalize raw lipidomics using MetaboAnalystR
# Follows EXACT R commands exported by MetaboAnalyst website:
#
#   InitDataObjects("pktable", "stat", FALSE, 150)
#   Read.TextData(..., "rowu", "disc")       # samples in rows
#   SanityCheckData
#   PerformSanityClosure
#   CheckContainsBlank
#   FilterVariable (no filtering — passthrough)
#   PreparePrenormData
#   Normalization("MedianNorm", "LogNorm", "ParetoNorm")
#
# Input: raw CSVs with samples in COLUMNS (our format from 00b)
#   → script transposes to samples in ROWS before passing to MetaboAnalystR
#
# Outputs → results/03_normalization/
# -------------------------------------------------------

source("R/00_load_packages.R")

# Initialize MetaboAnalystR globals before loading
.init_metabo_globals <- function() {
  globs <- list(
    default.dpi = 150, current.msg = "", msg.vec = character(0),
    err.vec = character(0), err.msg = "", norm.msg = "", norm.warn = "",
    sys.msg.vec = character(0), fig.count = 0, table.count = 0,
    imgName = "", filenm = "", rpath = "../../", anal.type = "",
    mSetObj = NULL, on.public.web = FALSE, input_filename = "",
    rawFileNms = character(0), rawfilenms.vec = character(0),
    rawClassNms = character(0), colVec = character(0),
    shapeVec = character(0), pca.cex = 1, pls.cex = 1,
    opls.cex = 1, spls.cex = 1, BHth = 0.05,
    module.count = 0, moduleNms.vec = character(0),
    fullUserPath = getwd()
  )
  for (nm in names(globs)) assign(nm, globs[[nm]], envir = .GlobalEnv)
}
.init_metabo_globals()
suppressMessages(suppressWarnings(library(MetaboAnalystR)))

# Patch InitDataObjects default.dpi self-reference bug
.patch_init <- function() {
  fn_body <- deparse(body(MetaboAnalystR::InitDataObjects))
  fn_body <- sub(
    'assign("default.dpi", default.dpi, envir = .GlobalEnv)',
    'assign("default.dpi", 150, envir = .GlobalEnv)',
    fn_body, fixed = TRUE
  )
  new_fn <- eval(parse(text = paste(
    c("function(data.type, anal.type, paired=FALSE, default.dpi=150) {",
      fn_body[-c(1, length(fn_body))], "}"),
    collapse = "\n"
  )))
  environment(new_fn) <- environment(MetaboAnalystR::InitDataObjects)
  utils::assignInNamespace("InitDataObjects", new_fn, ns = "MetaboAnalystR")
  assign("InitDataObjects", new_fn, envir = .GlobalEnv)
}
.patch_init()

cfg <- yaml::read_yaml("config/config.yml")
out <- cfg$output$s03
dir.create(out,                            showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$output$processed,           showWarnings = FALSE, recursive = TRUE)
dir.create("data/raw/lipids_normalized",   showWarnings = FALSE, recursive = TRUE)

# ── Load group assignments from trait file ────────────────
# Use Control/Intervention labels — consistent with trait data
trait_groups <- tryCatch({
  tr <- readRDS(file.path(cfg$output$processed, "traits_cleaned.rds"))
  setNames(
    ifelse(tr$group == 1, "Intervention", "Control"),
    rownames(tr)
  )
}, error = function(e) {
  message("  Note: trait file not found, using A/B dummy groups")
  NULL
})

# ── Helper: write temp CSV in "rowu" format ───────────────
# MetaboAnalyst website uses "rowu" = samples in rows
# Our raw files have samples in columns — we transpose here
write_rowu_csv <- function(raw_path, tmp_path, group_map = NULL) {
  raw <- read.csv(raw_path, header = FALSE, stringsAsFactors = FALSE,
                  check.names = FALSE)

  # Remove duplicate header row if present
  if (is.na(raw[1, 1]) || trimws(raw[1, 1]) == "") raw <- raw[-1, ]

  # Row 1 = Label (aliquot IDs)
  # Row 2 = #CLASS (participant IDs)
  # Rows 3+ = lipid name + intensities
  sample_ids      <- as.character(raw[1, -1])
  participant_ids <- as.character(raw[2, -1])
  feat_names      <- as.character(raw[3:nrow(raw), 1])

  int_mat <- apply(raw[3:nrow(raw), -1], 2,
                   function(x) as.numeric(as.character(x)))
  rownames(int_mat) <- feat_names
  colnames(int_mat) <- sample_ids

  # Transpose to samples x features (rowu format)
  int_mat_t <- t(int_mat)

  # Assign groups: Control/Intervention from trait file, or A/B fallback
  n_samp <- nrow(int_mat_t)
  if (!is.null(group_map)) {
    groups <- ifelse(participant_ids %in% names(group_map),
                     group_map[participant_ids],
                     "Control")
  } else {
    groups <- rep(c("Control", "Intervention"), length.out = n_samp)
  }

  # rowu format matching exactly what MetaboAnalyst website receives:
  # Row 1 (header): "POD GROWell ID", "ExperimentalGroup", feat1, feat2, ...
  # Row 2+ (data):  participantID,    Control/Intervention, val1,  val2, ...
  #
  # MetaboAnalystR reads:
  #   col 1 → sample names (participant IDs)
  #   col 2 → class labels (Control/Intervention)
  #   cols 3+ → feature intensities
  header <- c("POD GROWell ID", "ExperimentalGroup", feat_names)

  con <- file(tmp_path, "w")
  writeLines(paste(paste0('"', header, '"'), collapse = ","), con)
  for (i in seq_len(nrow(int_mat_t))) {
    row_vals <- c(participant_ids[i], groups[i], int_mat_t[i, ])
    writeLines(paste(paste0('"', row_vals, '"'), collapse = ","), con)
  }
  close(con)

  invisible(list(sample_ids      = sample_ids,
                 participant_ids = participant_ids,
                 feat_names      = feat_names,
                 groups          = groups))
}

# ── Helper: normalize one timepoint ──────────────────────
normalize_timepoint <- function(raw_path, out_path, tp_label, out_dir) {

  message("\n--- Normalizing: ", tp_label, " ---")

  # Step 1: write temp file in rowu format with dummy A/B groups
  tmp_path <- tempfile(fileext = ".csv")
  ids <- write_rowu_csv(raw_path, tmp_path, group_map = trait_groups)

  # Also read raw values for QC plots
  raw <- read.csv(raw_path, header = FALSE, stringsAsFactors = FALSE,
                  check.names = FALSE)
  if (is.na(raw[1, 1]) || trimws(raw[1, 1]) == "") raw <- raw[-1, ]
  feat_names <- as.character(raw[3:nrow(raw), 1])
  int_mat <- apply(raw[3:nrow(raw), -1], 2,
                   function(x) as.numeric(as.character(x)))
  rownames(int_mat) <- feat_names
  colnames(int_mat) <- ids$sample_ids

  message(sprintf("  Features: %d | Samples: %d",
                  nrow(int_mat), ncol(int_mat)))

  # Step 2: run exact MetaboAnalyst website pipeline
  if (exists("mSet", envir = .GlobalEnv)) rm("mSet", envir = .GlobalEnv)

  mSet <- InitDataObjects("pktable", "stat", FALSE, 150)
  mSet <- Read.TextData(mSet, tmp_path, "rowu", "disc")
  mSet <- SanityCheckData(mSet)

  # Skip PerformSanityClosure, CheckContainsBlank, FilterVariable
  # — these are website UI steps not needed for pure normalization
  mSet <- MetaboAnalystR:::PreparePrenormData(mSet)

  # Step 3: normalize — exact website settings
  mSet <- Normalization(mSet,
                        "MedianNorm", "LogNorm", "ParetoNorm",
                        ratio = FALSE, ratioNum = 20)

  unlink(tmp_path)

  # Move MetaboAnalystR intermediate files out of working directory
  metabo_temps <- c("data_orig_0.qs", "data_orig.qs", "preproc.orig.qs",
                    "preproc.qs", "prenorm.qs", "row_norm.qs",
                    "complete_norm.qs", "raw_dataview.csv", "data_proc.qs")
  for (f in metabo_temps) {
    if (file.exists(f)) {
      dest <- file.path(cfg$output$processed, paste0(tp_label, "_", f))
      file.rename(f, dest)
    }
  }

  # Step 4: extract normalized matrix
  # Website output is samples x features — transpose to features x samples
  norm_mat_t <- mSet$dataSet$norm   # samples x features
  norm_mat   <- t(norm_mat_t)       # features x samples

  message(sprintf("  Normalized range: [%.3f, %.3f]",
                  min(norm_mat, na.rm = TRUE), max(norm_mat, na.rm = TRUE)))

  # Step 5: QC plots
  # 1. Boxplot before/after
  pdf(file.path(out_dir, paste0(tp_label, "_01_sample_normalization.pdf")),
      width = 14, height = 5)
  par(mfrow = c(1, 2), mar = c(8, 4, 3, 1))
  boxplot(int_mat,
          main = paste0(tp_label, " — Before"),
          ylab = "Raw intensity", col = "#4393C3",
          las = 2, cex.axis = 0.55, outline = FALSE)
  boxplot(norm_mat,
          main = paste0(tp_label, " — After"),
          ylab = "Normalized intensity", col = "#D6604D",
          las = 2, cex.axis = 0.55, outline = FALSE)
  dev.off()

  # 2. Density plot
  pdf(file.path(out_dir, paste0(tp_label, "_02_density_plot.pdf")),
      width = 10, height = 4)
  par(mfrow = c(1, 2))
  d_raw  <- density(as.vector(int_mat),  na.rm = TRUE)
  d_norm <- density(as.vector(norm_mat), na.rm = TRUE)
  plot(d_raw,  main = paste0(tp_label, " — Before"),
       col = "#4393C3", lwd = 2, bty = "l")
  polygon(d_raw,  col = adjustcolor("#4393C3", 0.2), border = NA)
  plot(d_norm, main = paste0(tp_label, " — After"),
       col = "#D6604D", lwd = 2, bty = "l")
  polygon(d_norm, col = adjustcolor("#D6604D", 0.2), border = NA)
  dev.off()

  # 3. PCA scores plot
  pca_res <- prcomp(t(norm_mat), scale. = FALSE)
  pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
  scores  <- as.data.frame(pca_res$x[, 1:2])
  scores$Sample <- rownames(scores)

  p_pca <- ggplot(scores, aes(x = PC1, y = PC2, label = Sample)) +
    geom_point(size = 3, alpha = 0.8, color = "#2166AC") +
    labs(title = paste0(tp_label, " — PCA (post-normalization)"),
         x = paste0("PC1 (", pct_var[1], "%)"),
         y = paste0("PC2 (", pct_var[2], "%)")) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(file.path(out_dir, paste0(tp_label, "_03_pca_scores.pdf")),
         plot = p_pca, width = 7, height = 6)

  # Step 6: save normalized CSV
  # Restore original participant IDs in #CLASS
  pid_map      <- setNames(ids$participant_ids, ids$sample_ids)
  sample_order <- colnames(norm_mat)
  restored_cls <- pid_map[sample_order]

  out_df <- rbind(
    c("Label",  sample_order),
    c("#CLASS", restored_cls),
    cbind(rownames(norm_mat), norm_mat)
  )
  write.table(out_df, out_path,
              sep = ",", row.names = FALSE, col.names = FALSE, quote = TRUE)
  message("  Saved: ", out_path)

  data.frame(
    Timepoint  = tp_label,
    N_features = nrow(norm_mat),
    N_samples  = ncol(norm_mat),
    Norm_min   = round(min(norm_mat, na.rm = TRUE), 3),
    Norm_max   = round(max(norm_mat, na.rm = TRUE), 3),
    stringsAsFactors = FALSE
  )
}

# ── Process all three timepoints ──────────────────────────
timepoints <- list(
  list(label = "Baseline",
       raw   = cfg$data$lipids_raw$baseline,
       norm  = cfg$data$lipids$baseline),
  list(label = "Wk36_38",
       raw   = cfg$data$lipids_raw$wk36,
       norm  = cfg$data$lipids$wk36),
  list(label = "Postpartum",
       raw   = cfg$data$lipids_raw$postpartum,
       norm  = cfg$data$lipids$postpartum)
)

qc_summary <- lapply(timepoints, function(tp) {
  normalize_timepoint(
    raw_path = tp$raw,
    out_path = tp$norm,
    tp_label = tp$label,
    out_dir  = out
  )
}) %>% bind_rows()

write.csv(qc_summary,
          file.path(out, "normalization_qc_summary.csv"), row.names = FALSE)
write.xlsx(qc_summary,
           file.path(out, "normalization_qc_summary.xlsx"))

message("\nNormalization QC summary:")
print(qc_summary)
message("\nScript 03 complete → ", out)