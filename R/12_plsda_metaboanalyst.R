# R/12_plsda_metaboanalyst.R
# -------------------------------------------------------
# PLS-DA using MetaboAnalystR — mirrors EXACT MetaboAnalyst
# website workflow to reproduce published biplot results.
#
# Key design decisions:
#   - Data is ALREADY normalized (median + log10 + Pareto in script 03)
#     so Normalization() is called with NULL/NULL/NULL (no-op passthrough)
#   - Uses PLSR.Anal() which matches the MetaboAnalyst website PLS-DA
#   - Permutation testing via mSet$analSet$plsr for Q2/R2 reporting
#   - Model statistics (R2, Q2, perm p) extracted and saved to Excel

# Outputs → results/12_plsda_metaboanalyst/
# -------------------------------------------------------

source("R/00_load_packages.R")
source("R/00_figure_theme.R")

# ── MetaboAnalystR globals (required before loading) ─────
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

# Patch InitDataObjects default.dpi self-reference bug (same as script 03)
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

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Config and output directories ────────────────────────
cfg <- yaml::read_yaml("config/config.yml")

out <- file.path("results", "12_plsda_metaboanalyst")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

fig_dir <- file.path("results", "figures", "final")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load data ─────────────────────────────────────────────
# multiExpr_ID: normalized lipid expression per timepoint
# status_labels: binary outcome columns (apo_status, ppwr_status)
multiExpr_ID  <- readRDS(file.path(cfg$output$processed, "multiExpr_ID.rds"))
status_labels <- readRDS(file.path(cfg$output$processed, "lipid_status_labels.rds"))

tp_names   <- names(multiExpr_ID)
tp_display <- c("Baseline", "Wk36_38", "Postpartum")

# Outcomes to run — must match column names in status_labels
# These produce the biplots used in the paper
outcomes <- c("apo_status", "ppwr_status")

# ── Helper: write temp CSV in MetaboAnalyst rowu format ──
# MetaboAnalyst website expects:
#   Row 1 header: "Sample", "Label", feat1, feat2, ...
#   Row 2+:       sampleID, groupLabel, val1, val2, ...
# Normalization is NULL/NULL/NULL since data is pre-normalized

write_metabo_csv <- function(expr_mat, group_vec, tmp_path) {

  # Align samples
  common <- intersect(rownames(expr_mat), names(group_vec))
  if (length(common) == 0) stop("No overlapping samples")

  mat  <- expr_mat[common, , drop = FALSE]
  grps <- group_vec[common]

  # ── Remove samples with missing group labels ──────────
  # NA group values become "NA" strings which MetaboAnalystR
  # reads as a third group, causing a fatal error
  valid <- !is.na(grps) & nchar(as.character(grps)) > 0
  mat   <- mat[valid, , drop = FALSE]
  grps  <- grps[valid]

  if (length(unique(grps)) < 2) {
    stop("Fewer than 2 groups after removing NA labels")
  }

  message(sprintf("    Groups after NA removal: %s",
    paste(names(table(grps)), table(grps), sep="=", collapse=", ")))

  # Build output
  header <- c("Sample", "Label", colnames(mat))
  con <- file(tmp_path, "w")
  writeLines(paste(paste0('"', header, '"'), collapse = ","), con)
  for (i in seq_len(nrow(mat))) {
    row_vals <- c(rownames(mat)[i], as.character(grps[i]), mat[i, ])
    writeLines(paste(paste0('"', row_vals, '"'), collapse = ","), con)
  }
  close(con)

  invisible(rownames(mat))
}

# ── Helper: extract PLS-DA model statistics ───────────────
extract_plsda_stats <- function(mSet, tp_label, outcome_name,
                                n_samples, n_lipids) {

  plsr_obj <- mSet$analSet$plsr

  # R2 and Q2 per component from the pls object
  # pls::R2() and pls::RMSEP() work on the fitted model
  if (is.null(plsr_obj)) {
    warning("No PLSR object found for ", tp_label, " / ", outcome_name)
    return(NULL)
  }

  # R2Y (Y-variance explained) — from pls package
  r2y <- tryCatch({
    r <- pls::R2(plsr_obj, estimate = "CV")$val
    # r is array [responses, components+1] — take first response, skip intercept
    as.numeric(r[1, , 1])[-1]   # drop intercept (column 1)
  }, error = function(e) rep(NA_real_, 2))

  # Q2 (cross-validated) — from pls package
  q2 <- tryCatch({
    q <- pls::R2(plsr_obj, estimate = "CV")$val
    as.numeric(q[1, , 1])[-1]
  }, error = function(e) rep(NA_real_, 2))

  # MetaboAnalystR stores its own summary
  # mSet$analSet$plsr.cls.res contains permutation info if run
  perm_info <- mSet$analSet$plsr.cls.res

  # Permutation p-values — MetaboAnalystR stores these after PLSR.Anal
  # when permutation is run internally
  perm_pQ2  <- if (!is.null(perm_info) && "pQ2"  %in% names(perm_info)) {
    perm_info$pQ2
  } else NA_real_

  perm_pR2Y <- if (!is.null(perm_info) && "pR2Y" %in% names(perm_info)) {
    perm_info$pR2Y
  } else NA_real_

  # Cross-validation accuracy
  cv_acc <- if (!is.null(perm_info) && "acc" %in% names(perm_info)) {
    perm_info$acc
  } else NA_real_

  # VIP scores — stored in mSet$analSet$plsr.vip
  vip_df <- tryCatch({
    vip <- mSet$analSet$plsr.vip
    if (!is.null(vip)) {
      data.frame(
        Lipid     = names(vip),
        VIP_Comp1 = as.numeric(vip),
        row.names = NULL,
        stringsAsFactors = FALSE
      ) %>% dplyr::arrange(dplyr::desc(VIP_Comp1))
    } else NULL
  }, error = function(e) NULL)

  list(
    summary = data.frame(
      Timepoint   = tp_label,
      Outcome     = outcome_name,
      N_samples   = n_samples,
      N_lipids    = n_lipids,
      N_comp      = 2,
      Method      = "MetaboAnalystR PLSR",
      perm_pR2Y   = perm_pR2Y,
      perm_pQ2    = perm_pQ2,
      CV_accuracy = cv_acc,
      stringsAsFactors = FALSE
    ),
    vip_df = vip_df
  )
}

# ── Core function: run one PLS-DA via MetaboAnalystR ─────
run_plsda_metaboanalyst <- function(expr_mat, group_vec,
                                     tp_label, outcome_name,
                                     out_dir) {

  # Align
  common <- intersect(rownames(expr_mat), names(group_vec))
  if (length(common) < 10) {
    message("  Skipping ", tp_label, " / ", outcome_name,
            ": n = ", length(common), " < 10")
    return(NULL)
  }

  n_samples <- length(common)
  n_lipids  <- ncol(expr_mat)

  message(sprintf("  MetaboAnalyst PLS-DA: %s × %s  (n=%d, p=%d)",
                  tp_label, outcome_name, n_samples, n_lipids))

  # Write temp CSV
  tmp_csv <- tempfile(fileext = ".csv")
  write_metabo_csv(expr_mat, group_vec, tmp_csv)

  # Output prefix for MetaboAnalystR figure files
  pfx <- file.path(out_dir, paste0("plsda_MA_", tp_label, "_", outcome_name))

  # Set working directory context for MetaboAnalystR file outputs
  # MetaboAnalystR writes figures to working directory by default
  old_wd <- getwd()
  tmp_wd <- file.path(out_dir, paste0("tmp_", tp_label, "_", outcome_name))
  dir.create(tmp_wd, showWarnings = FALSE, recursive = TRUE)
  setwd(tmp_wd)

  tryCatch({

    # ── Steps 1–10: Initialize and load (no normalization) ──
    mSet <- InitDataObjects("pktable", "stat", FALSE, 600)
    mSet <- Read.TextData(mSet, tmp_csv, "rowu", "disc")
    mSet <- SanityCheckData(mSet)
    mSet <- MetaboAnalystR:::PreparePrenormData(mSet)

    # Normalization: NULL/NULL/NULL = no-op (data already normalized)
    mSet <- Normalization(mSet, "NULL", "NULL", "NULL",
                          ratio = FALSE, ratioNum = 20)

    # ── Step 11: PLS-DA (EXACT MetaboAnalyst website command) ──
    mSet <- PLSR.Anal(mSet, reg = TRUE)

    # ── Steps 12–18: Plots (EXACT MetaboAnalyst website commands) ──

    # Pair summary plot
    mSet <- PlotPLSPairSummary(mSet, "pls_pair_", "pdf", 72,
                               width = NA, 5)

    # 2D scores plot (comp 1 vs 2, 95% confidence ellipse)
    mSet <- PlotPLS2DScore(mSet, "pls_score2d_", "pdf", 72,
                           width = NA, 1, 2, 0.95, 0, 0, "na")

    # Loadings plot
    mSet <- PlotPLSLoading(mSet, "pls_loading_", "pdf", 72,
                           width = NA, 1, 2)

    # VIP importance plot (top 15, Comp 1)
    mSet <- PlotPLS.Imp(mSet, "pls_imp_", "pdf", 72,
                        width = NA, "vip", "Comp. 1", 15, FALSE)

    # Biplot (EXACT website command — comp 1 vs 2, top 10 loadings)
    mSet <- PlotPLSBiplot(mSet, "pls_biplot_", "pdf", 72,
                          width = NA, 1, 2, 10)

    # ── Copy outputs to named files in main out_dir ──────────
    plot_files <- list.files(tmp_wd, pattern = "\\.pdf$", full.names = TRUE)
    for (f in plot_files) {
      base <- basename(f)
      # Replace generic prefix with tp_label + outcome prefix
      new_name <- sub("^pls_", paste0("plsMA_", tp_label, "_",
                                      outcome_name, "_"), base)
      file.copy(f, file.path(out_dir, new_name), overwrite = TRUE)
    }

    # Also save high-res TIFF of biplot for submission
    biplot_pdf <- file.path(tmp_wd, "pls_biplot_.pdf")
    if (file.exists(biplot_pdf) &&
        requireNamespace("magick", quietly = TRUE)) {
      img <- magick::image_read_pdf(biplot_pdf, density = 600)
      magick::image_write(
        img,
        path = paste0(pfx, "_biplot.tiff"),
        format = "tiff",
        compression = "lzw"
      )
      message("  TIFF biplot saved: ", basename(paste0(pfx, "_biplot.tiff")))
    }

    # ── Extract statistics ────────────────────────────────────
    stats <- extract_plsda_stats(mSet, tp_label, outcome_name,
                                  n_samples, n_lipids)

    # ── Save Excel with model stats + VIP ────────────────────
    wb <- openxlsx::createWorkbook()

    if (!is.null(stats$summary)) {
      openxlsx::addWorksheet(wb, "Model_Summary")
      openxlsx::writeData(wb, "Model_Summary", stats$summary)
    }

    # Also save full PLSR model summary from MetaboAnalystR
    plsr_sum <- tryCatch({
      as.data.frame(summary(mSet$analSet$plsr)$explvar)
    }, error = function(e) NULL)

    if (!is.null(plsr_sum)) {
      openxlsx::addWorksheet(wb, "PLSR_Explained_Variance")
      openxlsx::writeData(wb, "PLSR_Explained_Variance", plsr_sum,
                          rowNames = TRUE)
    }

    if (!is.null(stats$vip_df)) {
      openxlsx::addWorksheet(wb, "VIP_scores")
      openxlsx::writeData(wb, "VIP_scores", stats$vip_df)
    }

    # Classification results table from MetaboAnalystR
    cls_res <- mSet$analSet$plsr.cls.res
    if (!is.null(cls_res)) {
      cls_df <- tryCatch(as.data.frame(cls_res), error = function(e) NULL)
      if (!is.null(cls_df)) {
        openxlsx::addWorksheet(wb, "Classification_Results")
        openxlsx::writeData(wb, "Classification_Results", cls_df,
                            rowNames = TRUE)
      }
    }

    openxlsx::saveWorkbook(wb, paste0(pfx, ".xlsx"), overwrite = TRUE)
    message("  Excel saved: ", basename(paste0(pfx, ".xlsx")))

    setwd(old_wd)
    unlink(tmp_csv)

    list(
      mSet        = mSet,
      stats       = stats,
      tp          = tp_label,
      outcome     = outcome_name,
      n_samples   = n_samples
    )

  }, error = function(e) {
    setwd(old_wd)
    unlink(tmp_csv)
    message("  ERROR in ", tp_label, " / ", outcome_name, ": ",
            conditionMessage(e))
    NULL
  })
}

# ── Main loop ─────────────────────────────────────────────
results <- list()

for (i in seq_along(tp_names)) {
  expr <- multiExpr_ID[[tp_names[i]]]$data   # samples x lipids
  tp   <- tp_display[i]

  for (outcome in outcomes) {

    if (!outcome %in% colnames(status_labels)) {
      message("  Outcome not found: ", outcome, " — skipping")
      next
    }

    grp_vec <- setNames(
      as.character(status_labels[[outcome]]),
      rownames(status_labels)
    )

    key <- paste0(tp, "_", outcome)
    results[[key]] <- run_plsda_metaboanalyst(
      expr_mat     = expr,
      group_vec    = grp_vec,
      tp_label     = tp,
      outcome_name = outcome,
      out_dir      = out
    )
  }
}

# ── Aggregate summary table ───────────────────────────────
summary_all <- lapply(names(results), function(key) {
  res <- results[[key]]
  if (is.null(res) || is.null(res$stats)) {
    return(data.frame(Key = key, Status = "Failed or Skipped",
                      stringsAsFactors = FALSE))
  }
  cbind(data.frame(Key = key, Status = "Complete",
                   stringsAsFactors = FALSE),
        res$stats$summary)
}) %>% dplyr::bind_rows()

wb_sum <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_sum, "PLS-DA_Summary")
openxlsx::writeData(wb_sum, "PLS-DA_Summary", summary_all)
openxlsx::saveWorkbook(wb_sum,
  file.path(out, "plsda_MA_summary_all.xlsx"),
  overwrite = TRUE)

message("\nPLS-DA MetaboAnalystR Summary:")
print(summary_all)
message("\nScript 12_plsda_metaboanalyst complete → ", out)
message("Biplots saved to: ", out)
message("To use in manuscript: copy biplot PDFs/TIFFs to supplementary figures")