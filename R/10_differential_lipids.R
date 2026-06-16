# R/10_differential_lipids.R
# -------------------------------------------------------
# Differential lipid analysis: APO and PPWR_E
#
# Statistics:
#   - Uses NORMALIZED/PREPARED lipid matrix from multiExpr_ID.rds
#   - Same test logic as prior working script:
#       Wilcoxon rank-sum or t-test controlled by config$differential$test
#   - Same BH FDR logic:
#       q_value_BH = p.adjust(p_value, method = config$differential$fdr_method)
#
# Fold change:
#   - Uses RAW lipid intensities from config$data$lipids_raw
#   - Calculates:
#       mean_raw_ctrl
#       mean_raw_case
#       fold_change_raw = mean_raw_case / mean_raw_ctrl
#       log2FC_raw      = log2(fold_change_raw)
#   - Also keeps normalized-scale mean difference:
#       mean_diff_scaled = mean_case_scaled - mean_ctrl_scaled
#
# Volcano plots:
#   - x-axis: log2FC_raw
#   - y-axis: -log10(p_value) from normalized-data test
#   - nominal p < cutoff: colored by direction
#   - FDR q < cutoff: black outline and labelled
#   - subtitle reports both nominal and FDR counts
#
# Outputs:
#   results/10_differential/<mode>/
#   results/figures/panels/fig4/fig4_A.pdf  (APO)
#   results/figures/panels/fig5/fig5_A.pdf  (PPWR)
# -------------------------------------------------------

source("R/00_load_packages.R")
source("R/00_figure_theme.R")

cfg <- yaml::read_yaml("config/config.yml")

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Settings ──────────────────────────────────────────────
diff_cfg    <- cfg$differential
outcomes    <- diff_cfg$outcomes
test_method <- tolower(diff_cfg$test %||% "wilcoxon")
apply_fdr   <- isTRUE(diff_cfg$apply_fdr)
fdr_method  <- diff_cfg$fdr_method %||% "BH"
fdr_cutoff  <- as.numeric(diff_cfg$fdr_cutoff %||% 0.10)
raw_cutoff  <- as.numeric(diff_cfg$raw_p_cutoff %||% 0.05)
top_n       <- as.integer(diff_cfg$top_n_labels %||% 12)

sig_cutoff <- if (apply_fdr) fdr_cutoff else raw_cutoff

# ── Output subfolder based on primary mode ────────────────
subfolder <- if (apply_fdr) {
  paste0("fdr_BH_q", gsub("\\.", "", sprintf("%.2f", fdr_cutoff)))
} else {
  paste0("raw_p", gsub("\\.", "", sprintf("%.2f", raw_cutoff)))
}

out <- file.path(cfg$output$s10, subfolder)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

message("Outcomes: ", paste(outcomes, collapse = ", "))
message("Test: ", test_method)
message("Primary mode: ",
        if (apply_fdr) paste0("FDR (", fdr_method, ", q<", fdr_cutoff, ")")
        else paste0("Raw p<", raw_cutoff))
message("Statistics source: normalized/prepared multiExpr_ID.rds")
message("Fold-change source: raw lipid intensities from config$data$lipids_raw")
message("Output: ", out)

# ── Timepoint labels and paths ────────────────────────────
tp_names <- c(
  "Identified_Baseline",
  "Identified_TP36_38weeks",
  "Identified_Postpartum"
)

tp_labels <- c("Baseline", "Wk36_38", "Postpartum")

tp_display <- c(
  "10\u201316 weeks",
  "36\u201338 weeks",
  "Postpartum"
)

raw_paths <- c(
  Baseline   = cfg$data$lipids_raw$baseline,
  Wk36_38    = cfg$data$lipids_raw$wk36,
  Postpartum = cfg$data$lipids_raw$postpartum
)

# ── Load normalized/prepared data for statistics ──────────
multiExpr_ID <- readRDS(file.path(cfg$output$processed, "multiExpr_ID.rds"))
demographic  <- readRDS(file.path(cfg$output$processed, "demographic_data_bmi.rds"))

# ── Helpers for reading raw lipid files ───────────────────
bad_sample_ids <- function(ids) {
  ids <- as.character(ids)

  is.null(ids) ||
    length(ids) == 0 ||
    any(is.na(ids)) ||
    any(ids == "") ||
    any(ids == "NA") ||
    all(grepl("^NA", ids))
}

clean_ids <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA
  x[x == "NA"] <- NA
  x
}

read_raw_lipid_file <- function(path, demographic_ids = NULL, timepoint_label = "timepoint") {

  if (is.null(path) || !file.exists(path)) {
    stop("Raw lipid file not found for ", timepoint_label, ": ", path)
  }

  message("\nReading RAW lipid file for ", timepoint_label, ": ", path)

  raw <- read.csv(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )

  message("  Raw CSV dimensions: ", paste(dim(raw), collapse = " x "))

  # Expected structure:
  # Row 1 = Label + aliquot/sample IDs
  # Row 2 = #CLASS + participant IDs
  # Rows 3+ = lipid names + raw intensity values
  label_ids       <- clean_ids(unlist(raw[1, -1]))
  participant_ids <- clean_ids(unlist(raw[2, -1]))

  lipid_names <- as.character(unlist(raw[3:nrow(raw), 1]))
  lipid_names <- make.unique(trimws(lipid_names))

  sample_ids <- participant_ids

  if (bad_sample_ids(sample_ids)) {
    message("  Row-2 participant IDs look invalid. Trying row-1 label IDs.")
    sample_ids <- label_ids
  }

  if (bad_sample_ids(sample_ids)) {
    if (!is.null(demographic_ids) && length(demographic_ids) == length(label_ids)) {
      message("  Both row-1 and row-2 IDs look invalid.")
      message("  Using demographic rownames because sample counts match.")
      sample_ids <- as.character(demographic_ids)
    } else {
      stop("Could not determine valid sample IDs for ", timepoint_label)
    }
  }

  sample_ids <- make.unique(as.character(sample_ids))

  int_mat <- raw[3:nrow(raw), -1, drop = FALSE]
  int_mat <- data.frame(
    lapply(int_mat, function(x) as.numeric(as.character(x))),
    check.names = FALSE
  )

  rownames(int_mat) <- lipid_names
  colnames(int_mat) <- sample_ids

  # Transpose to samples x lipids
  df_t <- as.data.frame(t(int_mat), check.names = FALSE)
  rownames(df_t) <- make.unique(trimws(as.character(rownames(df_t))))

  if (!is.null(demographic_ids)) {
    demographic_ids <- as.character(demographic_ids)
    overlap <- intersect(rownames(df_t), demographic_ids)

    message("  Overlap with demographic IDs: ", length(overlap))

    if (length(overlap) == 0) {
      if (nrow(df_t) == length(demographic_ids)) {
        message("  No ID overlap, but row counts match.")
        message("  Assigning demographic rownames to raw lipid data by order.")
        rownames(df_t) <- demographic_ids
      } else {
        stop(
          "No overlap between raw lipid sample IDs and demographic rownames for ",
          timepoint_label,
          ", and row counts do not match."
        )
      }
    } else {
      keep_ids <- demographic_ids[demographic_ids %in% rownames(df_t)]
      df_t <- df_t[keep_ids, , drop = FALSE]
    }
  }

  message("  RAW transposed/aligned dimensions: ", paste(dim(df_t), collapse = " x "))
  df_t
}

message("\nLoading raw lipid files for fold-change calculations...")
rawExpr_ID <- list()

for (i in seq_along(tp_labels)) {
  tp <- tp_labels[i]
  rawExpr_ID[[tp]] <- read_raw_lipid_file(
    path            = raw_paths[[tp]],
    demographic_ids = rownames(demographic),
    timepoint_label = tp
  )
}

# ── Safe fold-change helper ───────────────────────────────
calc_raw_fc <- function(raw_vals, grp, stat_valid) {

  # Use the same samples that were valid for the normalized-data statistical test.
  raw_vals <- raw_vals[stat_valid]
  grp      <- grp[stat_valid]

  raw_valid <- complete.cases(raw_vals, grp)
  raw_vals  <- as.numeric(raw_vals[raw_valid])
  grp       <- grp[raw_valid]

  g0_raw <- raw_vals[grp == 0]
  g1_raw <- raw_vals[grp == 1]

  mean_raw_ctrl <- mean(g0_raw, na.rm = TRUE)
  mean_raw_case <- mean(g1_raw, na.rm = TRUE)

  # True fold-change requires positive group means.
  if (is.na(mean_raw_ctrl) || is.na(mean_raw_case) ||
      mean_raw_ctrl <= 0 || mean_raw_case <= 0) {
    fold_change_raw <- NA_real_
    log2FC_raw      <- NA_real_
  } else {
    fold_change_raw <- mean_raw_case / mean_raw_ctrl
    log2FC_raw      <- log2(fold_change_raw)
  }

  list(
    mean_raw_ctrl   = mean_raw_ctrl,
    mean_raw_case   = mean_raw_case,
    fold_change_raw = fold_change_raw,
    log2FC_raw      = log2FC_raw,
    n_ctrl_raw_fc   = sum(grp == 0),
    n_case_raw_fc   = sum(grp == 1)
  )
}

# ── Test function: stats from normalized data, FC from raw ─
run_test <- function(expr_norm, expr_raw, group_vec, outcome_name, tp_label) {

  common <- Reduce(
    intersect,
    list(rownames(expr_norm), rownames(expr_raw), names(group_vec))
  )

  expr_norm <- expr_norm[common, , drop = FALSE]
  expr_raw  <- expr_raw[common, , drop = FALSE]
  group_vec <- group_vec[common]

  common_lipids <- intersect(colnames(expr_norm), colnames(expr_raw))

  missing_raw <- setdiff(colnames(expr_norm), colnames(expr_raw))
  if (length(missing_raw) > 0) {
    message("  ", tp_label, " x ", outcome_name,
            ": ", length(missing_raw),
            " normalized lipids missing from raw data; excluded from differential table.")
  }

  results <- lapply(common_lipids, function(lipid) {

    vals_norm <- as.numeric(expr_norm[[lipid]])
    vals_raw  <- as.numeric(expr_raw[[lipid]])
    grp       <- group_vec

    # Keep old statistics: validity based on normalized data + group only
    valid <- complete.cases(vals_norm, grp)
    vals_stat <- vals_norm[valid]
    grp_stat  <- grp[valid]

    if (length(unique(grp_stat)) < 2 || sum(valid) < 6) {
      return(NULL)
    }

    g0 <- vals_stat[grp_stat == 0]
    g1 <- vals_stat[grp_stat == 1]

    if (length(g0) < 2 || length(g1) < 2) {
      return(NULL)
    }

    tst <- if (test_method == "wilcoxon") {
      tryCatch(
        stats::wilcox.test(g0, g1),
        error = function(e) NULL
      )
    } else if (test_method %in% c("t", "ttest", "t-test")) {
      tryCatch(
        stats::t.test(g0, g1),
        error = function(e) NULL
      )
    } else {
      stop("Unsupported differential.test in config: ", test_method,
           ". Use 'wilcoxon' or 't-test'.")
    }

    if (is.null(tst)) {
      return(NULL)
    }

    fc <- calc_raw_fc(
      raw_vals   = vals_raw,
      grp        = grp,
      stat_valid = valid
    )

    mean_ctrl_scaled <- mean(g0, na.rm = TRUE)
    mean_case_scaled <- mean(g1, na.rm = TRUE)

    data.frame(
      Lipid       = lipid,
      Outcome     = outcome_name,
      Timepoint   = tp_label,

      # Statistical test result from normalized/prepared values
      p_value     = tst$p.value,

      # Raw-intensity effect size for volcano/reporting
      log2FC_raw      = fc$log2FC_raw,
      fold_change_raw = fc$fold_change_raw,
      mean_raw_ctrl   = fc$mean_raw_ctrl,
      mean_raw_case   = fc$mean_raw_case,

      # Normalized-scale summary retained for transparency
      mean_diff_scaled = mean_case_scaled - mean_ctrl_scaled,
      mean_ctrl_scaled = mean_ctrl_scaled,
      mean_case_scaled = mean_case_scaled,

      # Sample sizes
      n_ctrl        = sum(grp_stat == 0),
      n_case        = sum(grp_stat == 1),
      n_ctrl_raw_fc = fc$n_ctrl_raw_fc,
      n_case_raw_fc = fc$n_case_raw_fc,

      row.names = NULL
    )
  }) |>
    dplyr::bind_rows()

  if (is.null(results) || nrow(results) == 0) {
    return(NULL)
  }

  results[order(results$p_value), , drop = FALSE]
}

# ── Run differential tests ────────────────────────────────
message("\nRunning differential tests...")

diff_results_all <- list()

for (outcome in outcomes) {

  if (!(outcome %in% colnames(demographic))) {
    message("  Skipping '", outcome, "' — not found in demographic data")
    next
  }

  grp_vec <- setNames(demographic[[outcome]], rownames(demographic))
  grp_vec <- grp_vec[!is.na(grp_vec)]

  for (i in seq_along(tp_names)) {

    tp       <- tp_labels[i]
    expr_norm <- multiExpr_ID[[tp_names[i]]]$data
    expr_raw  <- rawExpr_ID[[tp]]

    res <- run_test(
      expr_norm    = expr_norm,
      expr_raw     = expr_raw,
      group_vec    = grp_vec,
      outcome_name = outcome,
      tp_label     = tp
    )

    if (is.null(res) || nrow(res) == 0) {
      message("  ", tp, " x ", outcome, ": no valid lipid tests")
      next
    }

    # Same FDR statistics as prior working script
    res$q_value_BH <- stats::p.adjust(res$p_value, method = fdr_method)

    # Explicit significance flags
    res$nominal_sig <- res$p_value < raw_cutoff
    res$fdr_sig     <- res$q_value_BH < fdr_cutoff

    # Legacy/significance column for downstream compatibility
    res$significant <- if (apply_fdr) res$fdr_sig else res$nominal_sig

    # Direction based on raw fold-change when available.
    # Fallback to scaled mean difference only if raw log2FC is NA.
    res$direction <- dplyr::case_when(
      !is.na(res$log2FC_raw) & res$log2FC_raw > 0 ~ "UP",
      !is.na(res$log2FC_raw) & res$log2FC_raw < 0 ~ "DOWN",
      is.na(res$log2FC_raw) & res$mean_diff_scaled > 0 ~ "UP",
      is.na(res$log2FC_raw) & res$mean_diff_scaled < 0 ~ "DOWN",
      TRUE ~ "NO_CHANGE"
    )

    key <- paste0(tp, "_", outcome)
    diff_results_all[[key]] <- res

    openxlsx::write.xlsx(
      res,
      file.path(out, paste0("differential_", key, ".xlsx")),
      overwrite = TRUE
    )

    message(sprintf(
      "  %s x %s: %d nominal (p<%.2f), %d FDR (q<%.2f) | total=%d | raw log2FC NA=%d",
      tp, outcome,
      sum(res$nominal_sig, na.rm = TRUE), raw_cutoff,
      sum(res$fdr_sig, na.rm = TRUE), fdr_cutoff,
      nrow(res),
      sum(is.na(res$log2FC_raw))
    ))
  }
}

# Save RDS with subfolder tag so downstream scripts know which version was used
saveRDS(
  diff_results_all,
  file.path(
    cfg$output$processed,
    paste0("differential_lipids_", subfolder, ".rds")
  )
)

# Also save as generic for downstream compatibility
saveRDS(
  diff_results_all,
  file.path(cfg$output$processed, "differential_lipids_all.rds")
)

# ── Volcano plots ─────────────────────────────────────────
make_volcano <- function(outcome, fig_id, panel_id, pos_label, neg_label) {

  keys <- paste0(tp_labels, "_", outcome)
  keys <- keys[keys %in% names(diff_results_all)]

  if (length(keys) == 0) {
    message("  No results for '", outcome, "' — skipping volcano")
    return(invisible(NULL))
  }

  plots <- lapply(seq_along(keys), function(i) {

    res  <- diff_results_all[[keys[i]]]
    tp_d <- tp_display[i]

    if (is.null(res) || nrow(res) == 0) {
      return(ggplot2::ggplot() + ggplot2::theme_void())
    }

    # Defensive checks in case objects are reused later
    if (!"q_value_BH" %in% colnames(res) || all(is.na(res$q_value_BH))) {
      res$q_value_BH <- stats::p.adjust(res$p_value, method = fdr_method)
    }
    if (!"nominal_sig" %in% colnames(res)) {
      res$nominal_sig <- res$p_value < raw_cutoff
    }
    if (!"fdr_sig" %in% colnames(res)) {
      res$fdr_sig <- res$q_value_BH < fdr_cutoff
    }
    if (!"direction" %in% colnames(res)) {
      res$direction <- ifelse(res$log2FC_raw > 0, "UP", "DOWN")
    }

    n_nominal <- sum(res$nominal_sig, na.rm = TRUE)
    n_fdr     <- sum(res$fdr_sig, na.rm = TRUE)

    # Label only FDR-significant lipids
    top_sig <- res |>
      dplyr::filter(fdr_sig) |>
      dplyr::arrange(q_value_BH, p_value) |>
      dplyr::slice_head(n = top_n)

    res$label <- ifelse(res$Lipid %in% top_sig$Lipid, res$Lipid, "")

    # Color nominally significant lipids by raw-FC direction
    res$color <- dplyr::case_when(
      res$nominal_sig & res$direction == "UP"   ~ "UP",
      res$nominal_sig & res$direction == "DOWN" ~ "DOWN",
      TRUE ~ "NS"
    )

    # For volcano x-axis, use raw log2FC. Lipids with NA raw log2FC
    # are omitted by ggplot from the x-axis plot, but remain in tables.
    p <- ggplot2::ggplot(
      res,
      ggplot2::aes(x = log2FC_raw, y = -log10(p_value))
    ) +

      # All lipids
      ggplot2::geom_point(
        ggplot2::aes(colour = color),
        size  = 0.9,
        alpha = 0.70,
        na.rm = TRUE
      ) +

      # FDR-significant lipids: black outline
      ggplot2::geom_point(
        data = res[res$fdr_sig & !is.na(res$log2FC_raw), , drop = FALSE],
        ggplot2::aes(x = log2FC_raw, y = -log10(p_value), fill = color),
        shape = 21,
        size  = 1.7,
        stroke = 0.35,
        colour = "black",
        alpha  = 0.95,
        show.legend = FALSE
      ) +

      ggplot2::scale_colour_manual(
        values = c(
          "UP"   = "#E63946",
          "DOWN" = "#1D3D8F",
          "NS"   = "#CCCCCC"
        ),
        labels = c(
          "UP"   = pos_label,
          "DOWN" = neg_label,
          "NS"   = "ns"
        ),
        name = NULL
      ) +

      ggplot2::scale_fill_manual(
        values = c(
          "UP"   = "#E63946",
          "DOWN" = "#1D3D8F",
          "NS"   = "#CCCCCC"
        ),
        guide = "none"
      ) +

      # Nominal p-value threshold line
      ggplot2::geom_hline(
        yintercept = -log10(raw_cutoff),
        linetype   = "dashed",
        colour     = GROWELL_COLORS$grey_dark,
        linewidth  = 0.3
      ) +

      ggplot2::labs(
        title = tp_d,
        subtitle = sprintf(
          "%d nominal p<%.2f; %d FDR q<%.2f",
          n_nominal, raw_cutoff, n_fdr, fdr_cutoff
        ),
        x = "Log2 fold change",
        y = "-log10(p)"
      ) +

      theme_growell(grid = "both") +

      ggplot2::theme(
        legend.position = "bottom",
        legend.text     = ggplot2::element_text(size = FS_SMALL),
        plot.title      = ggplot2::element_text(size = FS_BASE, face = "bold"),
        plot.subtitle   = ggplot2::element_text(
          size   = FS_SMALL,
          colour = GROWELL_COLORS$grey_dark
        )
      )

    if (requireNamespace("ggrepel", quietly = TRUE) && nrow(top_sig) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = res[res$label != "" & !is.na(res$log2FC_raw), , drop = FALSE],
        ggplot2::aes(label = label),
        size          = FS_SMALL / ggplot2::.pt,
        max.overlaps  = 20,
        segment.size  = 0.2,
        show.legend   = FALSE
      )
    }

    p
  })

  combined <- patchwork::wrap_plots(plots, nrow = 1) +
    patchwork::plot_annotation(
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size   = FS_TITLE,
          face   = "bold",
          family = FONT_FAMILY
        )
      )
    )

  ggplot2::ggsave(
    filename = file.path(out, paste0("volcano_", outcome, ".pdf")),
    plot     = combined,
    width    = FIG_WIDTH_FULL,
    height   = 80,
    units    = "mm",
    device   = grDevices::cairo_pdf
  )

  ggplot2::ggsave(
    filename = file.path(out, paste0("volcano_", outcome, ".png")),
    plot     = combined,
    width    = FIG_WIDTH_FULL,
    height   = 80,
    units    = "mm",
    dpi      = FIGURE_DPI
  )

  save_panel(
    combined,
    fig_id,
    panel_id,
    width_mm  = FIG_WIDTH_FULL,
    height_mm = 80
  )

  message("  Volcano: ", fig_id, "_", panel_id, " (", outcome, ")")
  invisible(combined)
}

message("\nGenerating volcano plots...")

make_volcano(
  outcome   = "apo",
  fig_id    = "fig4",
  panel_id  = "A",
  pos_label = "Higher in APO",
  neg_label = "Lower in APO"
)

make_volcano(
  outcome   = "ppwr_e",
  fig_id    = "fig5",
  panel_id  = "A",
  pos_label = "Higher in PPWR",
  neg_label = "Lower in PPWR"
)

# ── Summary table ─────────────────────────────────────────
summary_rows <- lapply(names(diff_results_all), function(key) {

  res <- diff_results_all[[key]]

  if (is.null(res) || nrow(res) == 0) {
    return(NULL)
  }

  if (!"nominal_sig" %in% colnames(res)) {
    res$nominal_sig <- res$p_value < raw_cutoff
  }
  if (!"fdr_sig" %in% colnames(res)) {
    res$fdr_sig <- res$q_value_BH < fdr_cutoff
  }

  data.frame(
    Analysis        = key,
    Total           = nrow(res),

    Sig_nominal_p05 = sum(res$nominal_sig, na.rm = TRUE),
    Sig_FDR_q10     = sum(res$fdr_sig,     na.rm = TRUE),

    Nominal_UP      = sum(res$nominal_sig & res$direction == "UP",
                          na.rm = TRUE),
    Nominal_DOWN    = sum(res$nominal_sig & res$direction == "DOWN",
                          na.rm = TRUE),

    FDR_UP          = sum(res$fdr_sig & res$direction == "UP",
                          na.rm = TRUE),
    FDR_DOWN        = sum(res$fdr_sig & res$direction == "DOWN",
                          na.rm = TRUE),

    Raw_log2FC_NA   = sum(is.na(res$log2FC_raw)),
    Raw_p_cutoff    = raw_cutoff,
    FDR_q_cutoff    = fdr_cutoff,
    FDR_method      = fdr_method,

    stringsAsFactors = FALSE
  )
}) |>
  dplyr::bind_rows()

openxlsx::write.xlsx(
  summary_rows,
  file.path(out, "differential_summary.xlsx"),
  overwrite = TRUE
)

print(summary_rows)

message("\nScript 10 complete → ", out)
message("Mode: ", subfolder)
message("Statistics were run on normalized/prepared lipid values.")
message("Fold-change columns and volcano x-axis use raw lipid intensities.")
message("To switch primary mode: change differential.apply_fdr in config/config.yml")

# # R/10_differential_lipids.R
# # -------------------------------------------------------
# # Differential lipid analysis: APO and PPWR_E
# #
# # Test:   Wilcoxon rank-sum or t-test (config differential.test)
# # FDR:    toggle via config differential.apply_fdr
# #   true  → BH correction, results → results/07_differential/fdr_BH_q{cutoff}/
# #   false → raw p-values,  results → results/07_differential/raw_p{cutoff}/
# #
# # Volcano plots:
# #   Figure 4A → results/figures/panels/fig4/fig4_A.pdf  (APO)
# #   Figure 5A → results/figures/panels/fig5/fig5_A.pdf  (PPWR)
# # -------------------------------------------------------
# source("R/00_load_packages.R")
# source("R/00_figure_theme.R")
# cfg <- yaml::read_yaml("config/config.yml")

# # ── Settings ──────────────────────────────────────────────
# diff_cfg    <- cfg$differential
# outcomes    <- diff_cfg$outcomes
# test_method <- diff_cfg$test
# apply_fdr   <- isTRUE(diff_cfg$apply_fdr)
# fdr_method  <- diff_cfg$fdr_method
# fdr_cutoff  <- diff_cfg$fdr_cutoff
# raw_cutoff  <- diff_cfg$raw_p_cutoff
# sig_cutoff  <- if (apply_fdr) fdr_cutoff else raw_cutoff
# top_n       <- diff_cfg$top_n_labels

# # ── Output subfolder based on method ─────────────────────
# subfolder <- if (apply_fdr) {
#   paste0("fdr_BH_q", gsub("\\.", "", sprintf("%.2f", fdr_cutoff)))
# } else {
#   paste0("raw_p", gsub("\\.", "", sprintf("%.2f", raw_cutoff)))
# }
# out <- file.path(cfg$output$s10, subfolder)
# dir.create(out, showWarnings = FALSE, recursive = TRUE)

# message("Outcomes: ", paste(outcomes, collapse = ", "))
# message("Test: ", test_method)
# message("Mode: ", if (apply_fdr) paste0("FDR (BH, q<", fdr_cutoff, ")") else paste0("Raw p<", raw_cutoff))
# message("Output: ", out)

# tp_names  <- c("Identified_Baseline","Identified_TP36_38weeks","Identified_Postpartum")
# tp_labels <- c("Baseline","Wk36_38","Postpartum")
# tp_display <- c("10\u201316 weeks","36\u201338 weeks","Postpartum")

# multiExpr_ID <- readRDS(file.path(cfg$output$processed, "multiExpr_ID.rds"))
# demographic  <- readRDS(file.path(cfg$output$processed, "demographic_data_bmi.rds"))

# # ── Test function ─────────────────────────────────────────
# run_test <- function(expr, group_vec, outcome_name, tp_label) {
#   common    <- intersect(rownames(expr), names(group_vec))
#   expr      <- expr[common, ]
#   group_vec <- group_vec[common]

#   results <- lapply(colnames(expr), function(lipid) {
#     vals  <- expr[[lipid]]
#     grp   <- group_vec
#     valid <- complete.cases(vals, grp)
#     vals  <- vals[valid]; grp <- grp[valid]
#     if (length(unique(grp)) < 2 || sum(valid) < 6) return(NULL)

#     g0 <- vals[grp == 0]; g1 <- vals[grp == 1]

#     tst <- if (test_method == "wilcoxon") {
#       tryCatch(wilcox.test(g0, g1), error = function(e) NULL)
#     } else {
#       tryCatch(t.test(g0, g1), error = function(e) NULL)
#     }
#     if (is.null(tst)) return(NULL)

#     data.frame(
#       Lipid     = lipid,
#       Outcome   = outcome_name,
#       Timepoint = tp_label,
#       log2FC    = mean(g1, na.rm = TRUE) - mean(g0, na.rm = TRUE),
#       p_value   = tst$p.value,
#       mean_ctrl = mean(g0, na.rm = TRUE),
#       mean_case = mean(g1, na.rm = TRUE),
#       n_ctrl    = sum(grp == 0),
#       n_case    = sum(grp == 1),
#       row.names = NULL
#     )
#   }) %>% bind_rows()

#   results[order(results$p_value), ]
# }

# # ── Run tests ─────────────────────────────────────────────
# message("\nRunning differential tests...")
# diff_results_all <- list()

# for (outcome in outcomes) {
#   if (!(outcome %in% colnames(demographic))) {
#     message("  Skipping '", outcome, "' — not found in demographic data")
#     next
#   }
#   grp_vec <- setNames(demographic[[outcome]], rownames(demographic))
#   grp_vec <- grp_vec[!is.na(grp_vec)]

#   for (i in seq_along(tp_names)) {
#     tp   <- tp_labels[i]
#     expr <- multiExpr_ID[[tp_names[i]]]$data
#     res  <- run_test(expr, grp_vec, outcome, tp)
#     if (is.null(res) || nrow(res) == 0) next

#     # Always calculate BH q-values so volcano plots can display
#     # both nominal and FDR-supported signals.
#     res$q_value_BH <- p.adjust(res$p_value, method = fdr_method)

#     # Explicit significance flags
#     res$nominal_sig <- res$p_value < raw_cutoff
#     res$fdr_sig     <- res$q_value_BH < fdr_cutoff

#     # Keep legacy column for downstream compatibility
#     res$significant <- if (apply_fdr) res$fdr_sig else res$nominal_sig

#     res$direction <- ifelse(res$log2FC > 0, "UP", "DOWN")
#     key <- paste0(tp, "_", outcome)
#     diff_results_all[[key]] <- res
#     write.xlsx(res, file.path(out, paste0("differential_", key, ".xlsx")))

#     message(sprintf(
#       "  %s x %s: %d nominal (p<%.2f), %d FDR (q<%.2f) | total=%d",
#       tp, outcome,
#       sum(res$nominal_sig, na.rm = TRUE), raw_cutoff,
#       sum(res$fdr_sig, na.rm = TRUE), fdr_cutoff,
#       nrow(res)
#     ))
#   }
# }

# # Save RDS with subfolder tag so 08 knows which version was used
# saveRDS(diff_results_all,
#         file.path(cfg$output$processed,
#                   paste0("differential_lipids_", subfolder, ".rds")))
# # Also save as generic for downstream compatibility
# saveRDS(diff_results_all,
#         file.path(cfg$output$processed, "differential_lipids_all.rds"))

# # ── Volcano plots ─────────────────────────────────────────
# make_volcano <- function(outcome, fig_id, panel_id, pos_label, neg_label) {
#   keys <- paste0(tp_labels, "_", outcome)
#   keys <- keys[keys %in% names(diff_results_all)]
#   if (length(keys) == 0) {
#     message("  No results for '", outcome, "' — skipping volcano"); return(invisible(NULL))
#   }

#   plots <- lapply(seq_along(keys), function(i) {
#     res  <- diff_results_all[[keys[i]]]
#     tp_d <- tp_display[i]
#     if (is.null(res) || nrow(res) == 0) return(ggplot() + theme_void())

#         # Ensure flags exist, even if loading older objects
#     if (!"q_value_BH" %in% colnames(res) || all(is.na(res$q_value_BH))) {
#       res$q_value_BH <- p.adjust(res$p_value, method = fdr_method)
#     }
#     if (!"nominal_sig" %in% colnames(res)) {
#       res$nominal_sig <- res$p_value < raw_cutoff
#     }
#     if (!"fdr_sig" %in% colnames(res)) {
#       res$fdr_sig <- res$q_value_BH < fdr_cutoff
#     }

#     n_nominal <- sum(res$nominal_sig, na.rm = TRUE)
#     n_fdr     <- sum(res$fdr_sig, na.rm = TRUE)

#     # Label only FDR-significant lipids
#     top_sig <- res %>%
#       filter(fdr_sig) %>%
#       arrange(q_value_BH, p_value) %>%
#       slice_head(n = top_n)

#     res$label <- ifelse(res$Lipid %in% top_sig$Lipid, res$Lipid, "")

#     # Color all nominally significant lipids by direction
#     res$color <- case_when(
#       res$nominal_sig & res$log2FC > 0 ~ "UP",
#       res$nominal_sig & res$log2FC < 0 ~ "DOWN",
#       TRUE ~ "NS"
#     )

#     p <- ggplot(res, aes(x = log2FC, y = -log10(p_value))) +

#       # All lipids: grey if not nominal, red/blue if nominal
#       geom_point(
#         aes(colour = color),
#         size = 0.9,
#         alpha = 0.7
#       ) +

#       # FDR-significant lipids: black outline
#       geom_point(
#         data = res[res$fdr_sig, , drop = FALSE],
#         aes(x = log2FC, y = -log10(p_value), fill = color),
#         shape = 21,
#         size = 1.7,
#         stroke = 0.35,
#         colour = "black",
#         alpha = 0.95,
#         show.legend = FALSE
#       ) +

#       scale_colour_manual(
#         values = c("UP" = "#E63946", "DOWN" = "#1D3D8F", "NS" = "#CCCCCC"),
#         labels = c("UP" = pos_label, "DOWN" = neg_label, "NS" = "ns"),
#         name   = NULL
#       ) +

#       scale_fill_manual(
#         values = c("UP" = "#E63946", "DOWN" = "#1D3D8F", "NS" = "#CCCCCC"),
#         guide = "none"
#       ) +

#       # Nominal p-value threshold line
#       geom_hline(
#         yintercept = -log10(raw_cutoff),
#         linetype = "dashed",
#         colour = GROWELL_COLORS$grey_dark,
#         linewidth = 0.3
#       ) +

#       labs(
#         title = tp_d,
#         subtitle = sprintf(
#           "%d nominal p<%.2f; %d FDR q<%.2f",
#           n_nominal, raw_cutoff, n_fdr, fdr_cutoff
#         ),
#         x = "log\u2082 fold change",
#         y = "-log\u2081\u2080(p)"
#       ) +

#       theme_growell(grid = "both") +

#       theme(
#         legend.position = "bottom",
#         legend.text     = element_text(size = FS_SMALL),
#         plot.title      = element_text(size = FS_BASE, face = "bold"),
#         plot.subtitle   = element_text(
#           size = FS_SMALL,
#           colour = GROWELL_COLORS$grey_dark
#         )
#       )
#     if (requireNamespace("ggrepel", quietly = TRUE) && nrow(top_sig) > 0) {
#       p <- p + ggrepel::geom_text_repel(
#         data = res[res$label != "", , drop = FALSE],
#         aes(label = label),
#         size = FS_SMALL / ggplot2::.pt,
#         max.overlaps = 20,
#         segment.size = 0.2,
#         show.legend = FALSE
#       )
#     }
#     p
#   })

#   combined <- wrap_plots(plots, nrow = 1) +
#     plot_annotation(
#       theme = theme(plot.title = element_text(size = FS_TITLE, face = "bold",
#                                                family = FONT_FAMILY))
#     )

#   ggsave(file.path(out, paste0("volcano_", outcome, ".pdf")),
#          plot = combined, width = FIG_WIDTH_FULL, height = 80, units = "mm",
#          device = cairo_pdf)
#   ggsave(file.path(out, paste0("volcano_", outcome, ".png")),
#          plot = combined, width = FIG_WIDTH_FULL, height = 80, units = "mm",
#          dpi = FIGURE_DPI)

#   save_panel(combined, fig_id, panel_id,
#              width_mm = FIG_WIDTH_FULL, height_mm = 80)

#   message("  Volcano: ", fig_id, "_", panel_id, " (", outcome, ")")
#   invisible(combined)
# }

# message("\nGenerating volcano plots...")
# make_volcano("apo",    "fig4", "A", "Higher in APO",  "Lower in APO")
# make_volcano("ppwr_e", "fig5", "A", "Higher in PPWR", "Lower in PPWR")

# # ── Summary table ─────────────────────────────────────────
# summary_rows <- lapply(names(diff_results_all), function(key) {
#   res <- diff_results_all[[key]]
#   if (is.null(res) || nrow(res) == 0) return(NULL)

#   if (!"nominal_sig" %in% colnames(res)) {
#     res$nominal_sig <- res$p_value < raw_cutoff
#   }
#   if (!"fdr_sig" %in% colnames(res)) {
#     res$fdr_sig <- res$q_value_BH < fdr_cutoff
#   }

#   data.frame(
#     Analysis        = key,
#     Total           = nrow(res),
#     Sig_nominal_p05 = sum(res$nominal_sig, na.rm = TRUE),
#     Sig_FDR_q10     = sum(res$fdr_sig, na.rm = TRUE),
#     Nominal_UP      = sum(res$nominal_sig & res$direction == "UP", na.rm = TRUE),
#     Nominal_DOWN    = sum(res$nominal_sig & res$direction == "DOWN", na.rm = TRUE),
#     FDR_UP          = sum(res$fdr_sig & res$direction == "UP", na.rm = TRUE),
#     FDR_DOWN        = sum(res$fdr_sig & res$direction == "DOWN", na.rm = TRUE)
#   )
# }) %>% bind_rows()

# write.xlsx(summary_rows, file.path(out, "differential_summary.xlsx"))
# print(summary_rows)

# message("\nScript 10 complete → ", out)
# message("Mode: ", subfolder)
# message("To switch: change differential.apply_fdr in config/config.yml")
