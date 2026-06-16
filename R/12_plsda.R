# R/12_plsda.R
# -------------------------------------------------------
# PLS-DA with ropls backend
#
# Provides all model statistics needed for manuscript:
#   R2X, R2Y per component (and cumulative)
#   Q2 cumulative (cross-validated)
#   Permutation p-values (pR2Y, pQ2)
#   CV accuracy
#   VIP scores
#
# Biplot uses publication-quality ggplot2 with Helvetica,
# 170mm width, 600 DPI — consistent with journal spec.
#
# Statistics reported in manuscript:
#   "PLS-DA was performed using the ropls package (v1.x).
#    Model validity was assessed by Q2 > 0 and permutation
#    p-value < 0.05 (999 permutations, 7-fold CV)."
#
# Outputs → results/12_plsda/
#   plsda_<tp>_<outcome>.xlsx        — per-model stats
#   plsda_<tp>_<outcome>_biplot.pdf  — publication figure
#   plsda_<tp>_<outcome>_biplot.tiff — 600 DPI
#   plsda_summary_all.xlsx           — all 6 models combined
#   plsda_manuscript_table.xlsx      — formatted for manuscript
# -------------------------------------------------------

source("R/00_load_packages.R")
source("R/00_figure_theme.R")

library(ggplot2)
library(dplyr)

`%||%` <- function(a, b) if (!is.null(a)) a else b

cfg <- yaml::read_yaml("config/config.yml")
out <- cfg$output$s12
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ── Load data ─────────────────────────────────────────────
multiExpr_ID  <- readRDS(file.path(cfg$output$processed, "multiExpr_ID.rds"))
status_labels <- readRDS(file.path(cfg$output$processed,
                                    "lipid_status_labels.rds"))

tp_names  <- names(multiExpr_ID)
tp_labels <- c("Baseline", "Wk36_38", "Postpartum")
tp_display <- c("10-16 weeks", "36-38 weeks", "Postpartum")

# ── Parameters from config ────────────────────────────────
cfg_outcomes  <- cfg$plsda$outcomes
all_outcomes  <- colnames(status_labels)
outcomes <- if (!is.null(cfg_outcomes) && length(cfg_outcomes) > 0)
  intersect(cfg_outcomes, all_outcomes) else all_outcomes

ncomp  <- cfg$plsda$ncomp   %||% 2L
permI  <- cfg$plsda$permI   %||% 999L
crossI <- cfg$plsda$crossI  %||% 7L
n_top  <- cfg$plsda$n_top_loadings %||% 10L

message("Running PLS-DA for: ", paste(outcomes, collapse = ", "))
message("Components: ", ncomp, " | CV folds: ", crossI,
        " | Permutations: ", permI)

# ── Colour palette ────────────────────────────────────────
GROUP_COLORS <- c(
  "APO_NO"   = "#1D3D8F", "APO_YES"  = "#E63946",
  "PPWR_NO"  = "#1D3D8F", "PPWR_YES" = "#E63946",
  "0"        = "#1D3D8F", "1"        = "#E63946"
)

# ── Helper: extract all stats from fitted ropls object ───
extract_ropls_stats <- function(fit, tp_label, outcome_name,
                                 n_samples, n_lipids) {
  mdf <- fit@modelDF
  sdf <- fit@summaryDF

  # Per-component values
  r2x_cum <- as.numeric(mdf[["R2X(cum)"]])
  r2y_cum <- as.numeric(mdf[["R2Y(cum)"]])
  q2_cum  <- as.numeric(mdf[["Q2(cum)"]])

  r2x_per <- diff(c(0, r2x_cum))
  r2y_per <- diff(c(0, r2y_cum))

  # Values at 2 components specifically
  nc <- min(2, length(r2x_per))

  # Permutation p-values
  perm_pR2Y <- tryCatch(
    as.numeric(sdf[nrow(sdf), "pR2Y"]), error = function(e) NA_real_)
  perm_pQ2  <- tryCatch(
    as.numeric(sdf[nrow(sdf), "pQ2"]),  error = function(e) NA_real_)

  # CV accuracy from modelDF if available
  cv_accuracy <- tryCatch({
    acc_col <- grep("Accuracy|accuracy", colnames(mdf),
                    value = TRUE, ignore.case = TRUE)[1]
    if (!is.na(acc_col)) as.numeric(mdf[nc, acc_col]) else NA_real_
  }, error = function(e) NA_real_)

  # Model validity flag
  q2_valid <- !is.na(q2_cum[nc]) && q2_cum[nc] > 0
  pq2_sig  <- !is.na(perm_pQ2) && perm_pQ2 < 0.05
  valid    <- q2_valid && pq2_sig

  # Print to console
  message(sprintf("    R2X: Comp1=%.3f, Comp2=%.3f (cum=%.3f)",
                  r2x_per[1], r2x_per[nc], r2x_cum[nc]))
  message(sprintf("    R2Y: Comp1=%.3f, Comp2=%.3f (cum=%.3f)",
                  r2y_per[1], r2y_per[nc], r2y_cum[nc]))
  message(sprintf("    Q2 (cum at %d comp): %.3f", nc, q2_cum[nc]))
  message(sprintf("    pR2Y=%.3f | pQ2=%.3f | Valid=%s",
                  perm_pR2Y %||% NA, perm_pQ2 %||% NA,
                  ifelse(valid, "YES", "NO")))

  # Summary row for manuscript table
  summary_row <- data.frame(
    Timepoint    = tp_label,
    Outcome      = outcome_name,
    N_samples    = as.integer(n_samples),
    N_lipids     = as.integer(n_lipids),
    N_comp       = as.integer(nc),
    CV_folds     = as.integer(crossI),
    Perm_n       = as.integer(permI),
    R2X_Comp1    = round(r2x_per[1],  3),
    R2X_Comp2    = round(r2x_per[nc], 3),
    R2X_cum      = round(r2x_cum[nc], 3),
    R2Y_Comp1    = round(r2y_per[1],  3),
    R2Y_Comp2    = round(r2y_per[nc], 3),
    R2Y_cum      = round(r2y_cum[nc], 3),
    Q2_cum       = round(q2_cum[nc],  3),
    perm_pR2Y    = round(perm_pR2Y,   3),
    perm_pQ2     = round(perm_pQ2,    3),
    CV_accuracy  = round(cv_accuracy, 3),
    Model_valid  = valid,
    stringsAsFactors = FALSE
  )

  # Full per-component table
  comp_df <- data.frame(
    Component    = seq_len(length(r2x_per)),
    R2X_per_comp = round(r2x_per, 4),
    R2X_cum      = round(r2x_cum, 4),
    R2Y_per_comp = round(r2y_per, 4),
    R2Y_cum      = round(r2y_cum, 4),
    Q2_cum       = round(q2_cum,  4),
    stringsAsFactors = FALSE
  )

  list(
    summary_row = summary_row,
    comp_df     = comp_df,
    mdf         = as.data.frame(mdf),
    sdf         = as.data.frame(sdf)
  )
}

# ── Helper: publication biplot ────────────────────────────
make_plsda_biplot <- function(fit, group_vec, tp_display_label,
                               outcome_name, n_top = 10) {

  sc <- tryCatch(as.data.frame(ropls::getScoreMN(fit)),
                  error = function(e) NULL)
  if (is.null(sc) || ncol(sc) < 2) return(NULL)

  colnames(sc)[1:2] <- c("Comp1", "Comp2")
  sc$Sample <- rownames(sc)
  sc$Group  <- as.character(group_vec[rownames(sc)])
  sc$Group[is.na(sc$Group)] <- "Unknown"

  ld <- tryCatch(as.data.frame(ropls::getLoadingMN(fit)),
                  error = function(e) NULL)
  if (is.null(ld)) return(NULL)
  colnames(ld)[1:2] <- c("Comp1", "Comp2")
  ld$Lipid <- rownames(ld)

  vip_vals  <- ropls::getVipVn(fit)
  top_names <- names(sort(vip_vals, decreasing = TRUE))[
    seq_len(min(n_top, nrow(ld)))]
  ld_top <- ld[ld$Lipid %in% top_names, ]

  # Axis labels with R2X
  r2x_cum <- fit@modelDF[["R2X(cum)"]]
  r2x_per <- diff(c(0, as.numeric(r2x_cum)))
  xlab <- sprintf("Component 1 (R2X=%.1f%%)", r2x_per[1] * 100)
  ylab <- sprintf("Component 2 (R2X=%.1f%%)",
                  if (length(r2x_per) >= 2) r2x_per[2] * 100 else 0)

  # Scale arrows
  score_range <- max(abs(c(range(sc$Comp1), range(sc$Comp2))),
                     na.rm = TRUE)
  load_range  <- max(abs(c(range(ld_top$Comp1), range(ld_top$Comp2))),
                     na.rm = TRUE)
  arrow_scale <- if (load_range > 0) 0.72 * score_range / load_range else 1
  ld_top <- ld_top %>%
    mutate(x_end = Comp1 * arrow_scale,
           y_end = Comp2 * arrow_scale)

  # Group colours
  grp_levels <- sort(unique(sc$Group[sc$Group != "Unknown"]))
  col_vals   <- sapply(grp_levels, function(g)
    GROUP_COLORS[g] %||% "#999999")

  # Outcome display label
  outcome_display <- dplyr::case_when(
    grepl("apo",  outcome_name, ignore.case = TRUE) ~ "APO composite",
    grepl("ppwr", outcome_name, ignore.case = TRUE) ~ "Excessive PPWR",
    TRUE ~ outcome_name
  )

  p <- ggplot2::ggplot() +

    # 95% confidence ellipses
    ggplot2::stat_ellipse(
      data = sc[sc$Group != "Unknown", ],
      ggplot2::aes(x = Comp1, y = Comp2, colour = Group, fill = Group),
      type = "norm", level = 0.95,
      geom = "polygon", alpha = 0.08,
      linewidth = 0.35, linetype = "dashed"
    ) +

    ggplot2::geom_hline(yintercept = 0, colour = "#DDDDDD",
                         linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "#DDDDDD",
                         linewidth = 0.3) +

    # Loading arrows
    ggplot2::geom_segment(
      data = ld_top,
      ggplot2::aes(x = 0, y = 0, xend = x_end, yend = y_end),
      arrow     = ggplot2::arrow(length = ggplot2::unit(0.12, "cm"),
                                  type = "closed"),
      colour    = "#555555", linewidth = 0.45, alpha = 0.8
    ) +

    # Score points
    ggplot2::geom_point(
      data = sc,
      ggplot2::aes(x = Comp1, y = Comp2, colour = Group),
      size = 2.0, alpha = 0.85
    ) +

    # Loading labels — geom_text to avoid viewport error in pipeline
    ggplot2::geom_text(
      data = ld_top,
      ggplot2::aes(x = x_end, y = y_end, label = Lipid),
      size          = 2.2, family = "Helvetica",
      fontface      = "bold", colour = "#333333",
      check_overlap = TRUE,
      nudge_y       = 0.05 * score_range
    ) +

    ggplot2::scale_colour_manual(values = col_vals, name = NULL) +
    ggplot2::scale_fill_manual(  values = col_vals, guide = "none") +

    ggplot2::labs(
      x        = xlab,
      y        = ylab,
      title    = tp_display_label,
      subtitle = outcome_display
    ) +

    ggplot2::theme_minimal(base_family = "Helvetica", base_size = 8) +
    ggplot2::theme(
      panel.border     = ggplot2::element_rect(
                           colour = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = ggplot2::element_line(
                           colour = "#EEEEEE", linewidth = 0.2),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(
                           size = 8, face = "bold", family = "Helvetica"),
      plot.subtitle    = ggplot2::element_text(
                           size = 7, colour = "#555555",
                           family = "Helvetica"),
      axis.text        = ggplot2::element_text(
                           size = 7, face = "bold",
                           colour = "black", family = "Helvetica"),
      axis.title       = ggplot2::element_text(
                           size = 7, face = "bold", family = "Helvetica"),
      legend.position  = "bottom",
      legend.text      = ggplot2::element_text(
                           size = 7, family = "Helvetica"),
      plot.margin      = ggplot2::margin(4, 4, 4, 4, unit = "mm")
    )

  p
}

# ── Core fitting function ─────────────────────────────────
run_plsda_ropls <- function(expr, group_vec, outcome_name,
                             tp_label, tp_display_label) {

  common <- intersect(rownames(expr), names(group_vec))
  if (length(common) == 0) {
    message("  Skipping ", tp_label, " / ", outcome_name,
            ": no overlapping samples")
    return(NULL)
  }

  X <- expr[common, , drop = FALSE]
  Y <- factor(as.character(group_vec[common]))

  # Remove constant/all-NA columns
  valid_cols <- apply(X, 2, function(x)
    !all(is.na(x)) && var(x, na.rm = TRUE) > 0)
  X <- X[, valid_cols, drop = FALSE]

  # Complete cases
  cc <- stats::complete.cases(X, as.integer(Y))
  X  <- X[cc, , drop = FALSE]
  Y  <- Y[cc]

  # Remove NA group labels
  valid_grp <- !is.na(Y) & nchar(as.character(Y)) > 0
  X <- X[valid_grp, , drop = FALSE]
  Y <- Y[valid_grp]

  if (length(unique(Y)) < 2 || nrow(X) < 10) {
    message("  Skipping ", tp_label, " / ", outcome_name,
            ": insufficient data after filtering (n=", nrow(X), ")")
    return(NULL)
  }

  n_samples <- nrow(X)
  n_lipids  <- ncol(X)

  message(sprintf("\n  PLS-DA [ropls]: %s x %s  (n=%d, p=%d)",
                  tp_label, outcome_name, n_samples, n_lipids))
  message(sprintf("    Groups: %s",
    paste(names(table(Y)), table(Y), sep = "=", collapse = ", ")))

  # Fit model
  fit <- tryCatch(
    ropls::opls(
      x         = X,
      y         = Y,
      predI     = ncomp,
      permI     = permI,
      crossvalI = crossI,
      info.txtC = "none",
      fig.pdfC  = "none"
    ),
    error = function(e) {
      message("  ropls error: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fit)) return(NULL)

  # Extract statistics
  stats <- extract_ropls_stats(fit, tp_label, outcome_name,
                                n_samples, n_lipids)

  # VIP scores
  vip_vals <- ropls::getVipVn(fit)
  vip_df   <- data.frame(
    Lipid     = names(vip_vals),
    VIP_Comp1 = as.numeric(vip_vals),
    row.names = NULL,
    stringsAsFactors = FALSE
  ) %>% dplyr::arrange(dplyr::desc(VIP_Comp1))

  # Loadings
  ld_all <- tryCatch({
    ld <- as.data.frame(ropls::getLoadingMN(fit))
    colnames(ld) <- paste0("Comp", seq_len(ncol(ld)))
    ld$Lipid <- rownames(ld)
    ld[, c("Lipid", setdiff(names(ld), "Lipid"))]
  }, error = function(e) NULL)

  # Scores
  scores_df <- tryCatch({
    sc <- as.data.frame(ropls::getScoreMN(fit))
    sc$Group <- as.character(group_vec[rownames(sc)])
    sc
  }, error = function(e) NULL)

  # Output prefix
  pfx <- file.path(out, paste0("plsda_", tp_label, "_", outcome_name))

  # Biplot — saved separately to avoid viewport errors
  bp <- make_plsda_biplot(
    fit               = fit,
    group_vec         = group_vec,
    tp_display_label  = tp_display_label,
    outcome_name      = outcome_name,
    n_top             = n_top
  )

  if (!is.null(bp)) {
    # PDF
    ggplot2::ggsave(
      filename = paste0(pfx, "_biplot.pdf"),
      plot     = bp,
      width    = 85, height = 85, units = "mm",
      dpi      = 600,
      device   = grDevices::cairo_pdf
    )
    # TIFF
    ggplot2::ggsave(
      filename = paste0(pfx, "_biplot.tiff"),
      plot     = bp,
      width    = 85, height = 85, units = "mm",
      dpi      = 600,
      device   = function(filename, ...) grDevices::tiff(
        filename, ..., compression = "lzw",
        type = "cairo", family = "Helvetica")
    )
    message("  Biplot saved: ", basename(pfx), "_biplot.pdf/.tiff")
  }

  # Permutation plot (base R from ropls)
  grDevices::pdf(paste0(pfx, "_permutation.pdf"), width = 7, height = 5)
  tryCatch(ropls::plot(fit, typeVc = "permutation"),
           error = function(e) NULL)
  grDevices::dev.off()

  # Score plot (base R from ropls)
  grDevices::pdf(paste0(pfx, "_scores.pdf"), width = 7, height = 6)
  tryCatch(ropls::plot(fit, typeVc = "x-score"),
           error = function(e) NULL)
  grDevices::dev.off()

  # ── Excel output ──────────────────────────────────────
  wb <- openxlsx::createWorkbook()

  # Sheet 1: Model summary — all stats needed for manuscript
  openxlsx::addWorksheet(wb, "Model_Summary")
  openxlsx::writeData(wb, "Model_Summary", stats$summary_row)

  # Sheet 2: Per-component table (R2X, R2Y, Q2 for each component)
  openxlsx::addWorksheet(wb, "Per_Component")
  openxlsx::writeData(wb, "Per_Component", stats$comp_df)

  # Sheet 3: Full ropls modelDF
  openxlsx::addWorksheet(wb, "ropls_modelDF")
  openxlsx::writeData(wb, "ropls_modelDF", stats$mdf)

  # Sheet 4: Full ropls summaryDF (permutation results)
  openxlsx::addWorksheet(wb, "ropls_summaryDF")
  openxlsx::writeData(wb, "ropls_summaryDF", stats$sdf)

  # Sheet 5: VIP scores all lipids
  openxlsx::addWorksheet(wb, "VIP_scores")
  openxlsx::writeData(wb, "VIP_scores", vip_df)

  # Sheet 6: Top 20 VIPs highlighted
  openxlsx::addWorksheet(wb, "VIP_top20")
  openxlsx::writeData(wb, "VIP_top20", head(vip_df, 20))

  # Sheet 7: Loadings
  if (!is.null(ld_all)) {
    openxlsx::addWorksheet(wb, "Loadings")
    openxlsx::writeData(wb, "Loadings", ld_all)
  }

  # Sheet 8: Scores with group labels
  if (!is.null(scores_df)) {
    openxlsx::addWorksheet(wb, "Scores")
    openxlsx::writeData(wb, "Scores", scores_df, rowNames = TRUE)
  }

  openxlsx::saveWorkbook(wb, paste0(pfx, ".xlsx"), overwrite = TRUE)
  message("  Excel saved: ", basename(paste0(pfx, ".xlsx")))

  list(
    model       = fit,
    stats       = stats,
    summary_row = stats$summary_row,
    vip_df      = vip_df,
    tp          = tp_label,
    outcome     = outcome_name,
    n_samples   = n_samples,
    n_lipids    = n_lipids
  )
}

# ── Main loop ─────────────────────────────────────────────
plsda_results <- list()

for (i in seq_along(tp_names)) {
  expr    <- multiExpr_ID[[tp_names[i]]]$data
  tp      <- tp_labels[i]
  tp_disp <- tp_display[i]

  for (outcome in outcomes) {
    grp_vec <- setNames(
      as.character(status_labels[[outcome]]),
      rownames(status_labels)
    )
    key <- paste0(tp, "_", outcome)
    plsda_results[[key]] <- run_plsda_ropls(
      expr             = expr,
      group_vec        = grp_vec,
      outcome_name     = outcome,
      tp_label         = tp,
      tp_display_label = tp_disp
    )
  }
}

saveRDS(plsda_results,
        file.path(cfg$output$processed, "plsda_results.rds"))

# ── Aggregate summary table ───────────────────────────────
summary_all <- lapply(names(plsda_results), function(key) {
  res <- plsda_results[[key]]
  if (is.null(res)) {
    return(data.frame(Key = key, Status = "Failed or Skipped",
                      stringsAsFactors = FALSE))
  }
  cbind(data.frame(Key = key, Status = "Complete",
                   stringsAsFactors = FALSE),
        res$summary_row)
}) %>% dplyr::bind_rows()

# ── Manuscript-ready table ────────────────────────────────
# Formatted with column names matching manuscript table conventions
manuscript_table <- summary_all %>%
  dplyr::filter(Status == "Complete") %>%
  dplyr::mutate(
    Outcome_display = dplyr::case_when(
      grepl("apo",  Outcome, ignore.case = TRUE) ~ "APO composite",
      grepl("ppwr", Outcome, ignore.case = TRUE) ~ "Excessive PPWR",
      TRUE ~ Outcome
    ),
    Timepoint_display = dplyr::case_when(
      Timepoint == "Baseline"   ~ "10-16 weeks",
      Timepoint == "Wk36_38"   ~ "36-38 weeks",
      Timepoint == "Postpartum" ~ "Postpartum",
      TRUE ~ Timepoint
    ),
    # Format for display
    R2X_display  = sprintf("%.3f / %.3f", R2X_Comp1, R2X_Comp2),
    R2Y_display  = sprintf("%.3f / %.3f", R2Y_Comp1, R2Y_Comp2),
    Q2_display   = sprintf("%.3f", Q2_cum),
    pQ2_display  = ifelse(perm_pQ2 < 0.001, "<0.001",
                           sprintf("%.3f", perm_pQ2)),
    Valid_display = ifelse(Model_valid, "Yes", "No")
  ) %>%
  dplyr::select(
    Timepoint_display, Outcome_display,
    N_samples, R2X_display, R2Y_display,
    Q2_display, pQ2_display, Valid_display
  )

colnames(manuscript_table) <- c(
  "Timepoint", "Outcome", "N",
  "R2X (Comp1/Comp2)", "R2Y (Comp1/Comp2)",
  "Q2 (cumulative)", "p(Q2)", "Model valid"
)

# Save both tables
wb_sum <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb_sum, "Full_Summary")
openxlsx::writeData(wb_sum, "Full_Summary", summary_all)

openxlsx::addWorksheet(wb_sum, "Manuscript_Table")
openxlsx::writeData(wb_sum, "Manuscript_Table", manuscript_table)

# Style the manuscript table
hdr_style <- openxlsx::createStyle(
  fontName = "Helvetica", fontSize = 10,
  textDecoration = "bold",
  fgFill = "#2C3E50", fontColour = "#FFFFFF",
  halign = "CENTER", wrapText = TRUE)
body_style <- openxlsx::createStyle(
  fontName = "Helvetica", fontSize = 10)
valid_style <- openxlsx::createStyle(
  fontName = "Helvetica", fontSize = 10,
  fgFill = "#C8E6C9", fontColour = "#1B5E20",
  textDecoration = "bold")

openxlsx::addStyle(wb_sum, "Manuscript_Table", hdr_style,
                    rows = 1, cols = 1:ncol(manuscript_table))
openxlsx::addStyle(wb_sum, "Manuscript_Table", body_style,
                    rows = 2:(nrow(manuscript_table) + 1),
                    cols = 1:ncol(manuscript_table),
                    gridExpand = TRUE)

# Highlight valid models in green
valid_rows <- which(manuscript_table[["Model valid"]] == "Yes") + 1
if (length(valid_rows) > 0) {
  openxlsx::addStyle(wb_sum, "Manuscript_Table", valid_style,
                      rows = valid_rows,
                      cols = ncol(manuscript_table))
}

openxlsx::setColWidths(wb_sum, "Manuscript_Table",
                        cols = 1:ncol(manuscript_table),
                        widths = c(14, 16, 6, 18, 18, 16, 10, 12))

openxlsx::saveWorkbook(wb_sum,
  file.path(out, "plsda_summary_all.xlsx"),
  overwrite = TRUE)

message("\nPLS-DA ropls Summary:")
print(summary_all[, c("Key", "Status", "R2X_Comp1", "R2X_Comp2",
                       "R2Y_Comp1", "R2Y_Comp2", "Q2_cum",
                       "perm_pQ2", "Model_valid")])

message("\nManuscript Table:")
print(manuscript_table)

message("\nScript 12_plsda complete -> ", out)
message("Key outputs:")
message("  plsda_summary_all.xlsx — Full_Summary + Manuscript_Table sheets")
message("  plsda_<tp>_<outcome>.xlsx — per-model stats (8 sheets each)")
message("  plsda_<tp>_<outcome>_biplot.pdf/.tiff — publication figures")
# # R/12_plsda.R
# # -------------------------------------------------------
# # PLS-DA with ropls backend
# #   - Q², R²X / R²Y, permutation p-value (native ropls)
# #   - Biplot (scores + loadings overlay, ggplot2 + ggrepel)
# #   - Config-driven outcome selection via plsda.outcomes
# #   - Per-timepoint × outcome Excel summary
# # Outputs → results/09_plsda/
# # -------------------------------------------------------

# source("R/00_load_packages.R")
# source("R/00_figure_theme.R")

# cfg        <- yaml::read_yaml("config/config.yml")
# out        <- cfg$output$s12
# dir.create(out, showWarnings = FALSE, recursive = TRUE)

# multiExpr_ID  <- readRDS(file.path(cfg$output$processed, "multiExpr_ID.rds"))
# status_labels <- readRDS(file.path(cfg$output$processed, "lipid_status_labels.rds"))

# tp_names  <- names(multiExpr_ID)
# tp_labels <- c("Baseline", "Wk36_38", "Postpartum")

# # ── Outcomes: pull from config, fall back to all available ──────────────────
# cfg_outcomes <- cfg$plsda$outcomes          # character vector or NULL
# all_outcomes <- colnames(status_labels)
# outcomes <- if (!is.null(cfg_outcomes) && length(cfg_outcomes) > 0) {
#   missing_oc <- setdiff(cfg_outcomes, all_outcomes)
#   if (length(missing_oc) > 0)
#     warning("plsda.outcomes not found in status_labels: ",
#             paste(missing_oc, collapse = ", "))
#   intersect(cfg_outcomes, all_outcomes)
# } else {
#   all_outcomes
# }
# message("Running PLS-DA for outcomes: ", paste(outcomes, collapse = ", "))

# # ── PLS-DA parameters from config ───────────────────────────────────────────
# ncomp   <- cfg$plsda$ncomp      # number of components
# permI   <- cfg$plsda$permI      # permutation iterations (e.g. 999)
# crossI  <- cfg$plsda$crossI     # cross-validation folds (e.g. 7)
# n_top   <- cfg$plsda$n_top_loadings %||% 10L   # top loadings to show in biplot

# # ── Colour palette (group labels come from factor levels) ───────────────────
# group_colours <- cfg$plsda$group_colours   # named list from config, or NULL
# default_pal   <- c(
#   "#4DAFCA", "#E94D36", "#4DAF4A", "#984EA3",
#   "#FF7F00", "#A65628", "#F781BF", "#999999"
# )

# # ── Helper: %||% (null-coalescing) ───────────────────────────────────────────
# `%||%` <- function(a, b) if (!is.null(a)) a else b

# # ── Biplot function ──────────────────────────────────────────────────────────
# #' @param fit        ropls opls object (2-component PLS-DA)
# #' @param group_vec  named factor/character; names = sample IDs
# #' @param tp_label   string for plot title
# #' @param outcome_label string for plot subtitle
# #' @param n_top      number of top-VIP loadings to display
# #' @param pal        named character vector of colours keyed to group levels,
# #'                   or NULL to use defaults
# make_plsda_biplot <- function(fit, group_vec, tp_label, outcome_label,
#                                n_top = 10, pal = NULL) {

#   # --- scores -----------------------------------------------------------------
#   sc <- as.data.frame(ropls::getScoreMN(fit))
#   if (ncol(sc) < 2) {
#     message("  Biplot skipped: fewer than 2 components in model")
#     return(NULL)
#   }
#   colnames(sc)[1:2] <- c("Comp1", "Comp2")
#   sc$Sample <- rownames(sc)
#   sc$Group  <- as.character(group_vec[rownames(sc)])

#   # --- loadings (X-side) ------------------------------------------------------
#   ld <- as.data.frame(ropls::getLoadingMN(fit))
#   colnames(ld)[1:2] <- c("Comp1", "Comp2")
#   ld$Lipid <- rownames(ld)

#   # top n by VIP on Comp1
#   vip_vals  <- ropls::getVipVn(fit)
#   top_names <- names(sort(vip_vals, decreasing = TRUE))[seq_len(min(n_top, nrow(ld)))]
#   ld_top    <- ld[ld$Lipid %in% top_names, ]

#   # --- axis labels (R²X per component, matching MetaboAnalyst convention) ----
#   summ <- ropls::getSummaryDF(fit)
#   # modelDF rows correspond to each component; R2X(cum) is cumulative
#   # Back-calculate per-component R²X
#   r2x_cum <- fit@modelDF[["R2X(cum)"]]
#   r2x_per <- diff(c(0, r2x_cum))
#   xlab <- sprintf("Component 1 (%.1f%%)", r2x_per[1] * 100)
#   ylab <- if (length(r2x_per) >= 2)
#     sprintf("Component 2 (%.1f%%)", r2x_per[2] * 100)
#   else
#     "Component 2"

#   # --- scale arrows into score space ------------------------------------------
#   score_range <- max(abs(c(range(sc$Comp1), range(sc$Comp2))), na.rm = TRUE)
#   load_range  <- max(abs(c(range(ld_top$Comp1), range(ld_top$Comp2))), na.rm = TRUE)
#   arrow_scale <- if (load_range > 0) 0.75 * score_range / load_range else 1
#   ld_top <- ld_top %>%
#     dplyr::mutate(x_end = Comp1 * arrow_scale,
#                   y_end = Comp2 * arrow_scale)

#   # --- colour palette ---------------------------------------------------------
#   grp_levels <- sort(unique(sc$Group))
#   if (!is.null(pal) && all(grp_levels %in% names(pal))) {
#     col_vals <- pal[grp_levels]
#   } else {
#     col_vals <- setNames(default_pal[seq_along(grp_levels)], grp_levels)
#   }

#   # --- plot -------------------------------------------------------------------
#   ggplot2::ggplot() +
#     # sample scores
#     ggplot2::geom_point(
#       data = sc,
#       ggplot2::aes(x = Comp1, y = Comp2, colour = Group),
#       size = 2.5, alpha = 0.85
#     ) +
#     # loading arrows
#     ggplot2::geom_segment(
#       data = ld_top,
#       ggplot2::aes(x = 0, y = 0, xend = x_end, yend = y_end),
#       arrow     = ggplot2::arrow(length = ggplot2::unit(0.18, "cm"),
#                                   type = "closed"),
#       colour    = "grey30",
#       linewidth = 0.4
#     ) +
#     # lipid labels — ggrepel avoids overlap
#     ggrepel::geom_text_repel(
#       data = ld_top,
#       ggplot2::aes(x = x_end, y = y_end, label = Lipid),
#       size         = 2.8,
#       colour       = "grey15",
#       max.overlaps = 25,
#       segment.size = 0.2,
#       box.padding  = 0.35
#     ) +
#     ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
#                          colour = "grey70", linewidth = 0.3) +
#     ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
#                          colour = "grey70", linewidth = 0.3) +
#     ggplot2::scale_colour_manual(values = col_vals, name = NULL) +
#     ggplot2::labs(
#       x        = xlab,
#       y        = ylab,
#       title    = paste("PLS-DA Biplot:", tp_label),
#       subtitle = outcome_label
#     ) +
#     ggplot2::theme_bw(base_size = 11) +
#     ggplot2::theme(
#       legend.position  = "bottom",
#       panel.grid.minor = ggplot2::element_blank(),
#       plot.title       = ggplot2::element_text(face = "bold", size = 11),
#       plot.subtitle    = ggplot2::element_text(size = 9, colour = "grey40")
#     )
# }

# # ── Core fitting function ────────────────────────────────────────────────────
# #' Fit a PLS-DA model with ropls, export plots + Excel, return summary row
# #'
# #' @return list with slots: model, summary_row, vip_df, or NULL on failure
# run_plsda_ropls <- function(expr, group_vec, outcome_name, tp_label) {

#   # --- align samples ----------------------------------------------------------
#   common <- intersect(rownames(expr), names(group_vec))
#   if (length(common) == 0) {
#     message("  Skipping ", tp_label, " / ", outcome_name, ": no overlapping samples")
#     return(NULL)
#   }
#   X <- expr[common, , drop = FALSE]
#   Y <- factor(group_vec[common])

#   # --- filter constant / all-NA columns ---------------------------------------
#   valid_cols <- apply(X, 2, function(x) !all(is.na(x)) && var(x, na.rm = TRUE) > 0)
#   X <- X[, valid_cols, drop = FALSE]

#   # --- complete cases ---------------------------------------------------------
#   cc <- stats::complete.cases(X, as.integer(Y))
#   X <- X[cc, , drop = FALSE]
#   Y <- Y[cc]

#   if (length(unique(Y)) < 2) {
#     message("  Skipping ", tp_label, " / ", outcome_name,
#             ": outcome has < 2 levels after filtering")
#     return(NULL)
#   }
#   if (nrow(X) < 10) {
#     message("  Skipping ", tp_label, " / ", outcome_name,
#             ": n = ", nrow(X), " < 10 after filtering")
#     return(NULL)
#   }

#   message(sprintf("  PLS-DA [ropls]: %s × %s  (n=%d, p=%d, permI=%d)",
#                   tp_label, outcome_name, nrow(X), ncol(X), permI))

#   # --- fit ropls opls ---------------------------------------------------------
#   fit <- tryCatch(
#   ropls::opls(
#     x     = X,
#     y     = Y,
#     predI = ncomp,
#     permI = permI,
#     crossvalI = crossI,        # was crossI in older versions
#     info.txtC = "none",        # suppresses printed output
#     fig.pdfC  = "none"         # suppresses auto PDF generation
#   ),
#   error = function(e) {
#     message("  ropls error: ", conditionMessage(e))
#     NULL
#   }
# )
#   if (is.null(fit)) return(NULL)

#   # --- extract summary stats --------------------------------------------------
#   mdf      <- fit@modelDF                     # one row per component
#   r2x_cum  <- mdf[["R2X(cum)"]]
#   r2y_cum  <- mdf[["R2Y(cum)"]]
#   q2_cum   <- mdf[["Q2(cum)"]]
#   r2x_per  <- diff(c(0, r2x_cum))
#   r2y_per  <- diff(c(0, r2y_cum))

#   # permutation p-values live in summaryDF (one row per perm stat)
#   sdf      <- fit@summaryDF
#   perm_pR2Y <- if ("pR2Y" %in% colnames(sdf)) sdf[nrow(sdf), "pR2Y"] else NA_real_
#   perm_pQ2  <- if ("pQ2"  %in% colnames(sdf)) sdf[nrow(sdf), "pQ2"]  else NA_real_

#   # VIP scores (single vector, Comp1-focused from ropls)
#   vip_vals <- ropls::getVipVn(fit)
#   vip_df   <- data.frame(
#     Lipid     = names(vip_vals),
#     VIP_Comp1 = vip_vals,
#     row.names = NULL
#   ) %>% dplyr::arrange(dplyr::desc(VIP_Comp1))

#   # --- output prefix ----------------------------------------------------------
#   pfx <- file.path(out, paste0("plsda_", tp_label, "_", outcome_name))

#   # --- biplot (publication quality TIFF + PDF) --------------------------------
#   bp <- make_plsda_biplot(
#     fit            = fit,
#     group_vec      = group_vec,
#     tp_label       = tp_label,
#     outcome_label  = outcome_name,
#     n_top          = n_top,
#     pal            = group_colours
#   )
#   if (!is.null(bp)) {
#     ggplot2::ggsave(paste0(pfx, "_biplot.pdf"),  bp, width = 7,   height = 6,   dpi = 300)
#     ggplot2::ggsave(paste0(pfx, "_biplot.tiff"), bp, width = 7,   height = 6,   dpi = 600,
#                     compression = "lzw")
#   }

#   # --- permutation plot (ropls built-in, saved to PDF) ------------------------
#   pdf(paste0(pfx, "_permutation.pdf"), width = 7, height = 5)
#   tryCatch(ropls::plot(fit, typeVc = "permutation"), error = function(e) NULL)
#   dev.off()

#   # --- score plot (ropls built-in) --------------------------------------------
#   pdf(paste0(pfx, "_scores.pdf"), width = 7, height = 6)
#   tryCatch(ropls::plot(fit, typeVc = "x-score"), error = function(e) NULL)
#   dev.off()

#   # --- Excel output -----------------------------------------------------------
#   wb <- openxlsx::createWorkbook()

#   # Sheet 1: model summary
#   n_comp_fit <- nrow(mdf)
#   summary_sheet <- data.frame(
#     Timepoint   = tp_label,
#     Outcome     = outcome_name,
#     N_samples   = nrow(X),
#     N_lipids    = ncol(X),
#     N_comp      = n_comp_fit,
#     CrossI      = crossI,
#     PermI       = permI
#   )
#   for (k in seq_len(n_comp_fit)) {
#     summary_sheet[[paste0("R2X_Comp",  k)]] <- round(r2x_per[k],  4)
#     summary_sheet[[paste0("R2Y_Comp",  k)]] <- round(r2y_per[k],  4)
#     summary_sheet[[paste0("Q2_cum_Comp", k)]] <- round(q2_cum[k], 4)
#   }
#   summary_sheet$perm_pR2Y <- perm_pR2Y
#   summary_sheet$perm_pQ2  <- perm_pQ2

#   openxlsx::addWorksheet(wb, "Model_Summary")
#   openxlsx::writeData(wb, "Model_Summary", summary_sheet)

#   # Sheet 2: full modelDF from ropls
#   openxlsx::addWorksheet(wb, "modelDF")
#   openxlsx::writeData(wb, "modelDF", as.data.frame(mdf))

#   # Sheet 3: VIP scores
#   openxlsx::addWorksheet(wb, "VIP_scores")
#   openxlsx::writeData(wb, "VIP_scores", vip_df)

#   # Sheet 4: loadings
#   ld_all <- as.data.frame(ropls::getLoadingMN(fit))
#   colnames(ld_all) <- paste0("Comp", seq_len(ncol(ld_all)))
#   ld_all$Lipid <- rownames(ld_all)
#   ld_all <- ld_all[, c("Lipid", setdiff(names(ld_all), "Lipid"))]
#   openxlsx::addWorksheet(wb, "Loadings")
#   openxlsx::writeData(wb, "Loadings", ld_all)

#   openxlsx::saveWorkbook(wb, paste0(pfx, ".xlsx"), overwrite = TRUE)

#   # --- return -----------------------------------------------------------------
#   list(
#     model       = fit,
#     vip_df      = vip_df,
#     summary_row = summary_sheet,
#     tp          = tp_label,
#     outcome     = outcome_name,
#     n_samples   = nrow(X),
#     n_lipids    = ncol(X)
#   )
# }

# # ── Main loop ────────────────────────────────────────────────────────────────
# plsda_results <- list()

# for (i in seq_along(tp_names)) {
#   expr <- multiExpr_ID[[tp_names[i]]]$data
#   tp   <- tp_labels[i]

#   for (outcome in outcomes) {
#     key     <- paste0(tp, "_", outcome)
#     grp_vec <- setNames(status_labels[[outcome]], rownames(status_labels))

#     plsda_results[[key]] <- run_plsda_ropls(
#       expr         = expr,
#       group_vec    = grp_vec,
#       outcome_name = outcome,
#       tp_label     = tp
#     )
#   }
# }

# saveRDS(plsda_results, file.path(cfg$output$processed, "plsda_results.rds"))

# # ── Aggregate summary table ──────────────────────────────────────────────────
# summary_all <- lapply(names(plsda_results), function(key) {
#   res <- plsda_results[[key]]
#   if (is.null(res)) {
#     return(data.frame(Key = key, Status = "Failed or Skipped",
#                       stringsAsFactors = FALSE))
#   }
#   cbind(data.frame(Key = key, Status = "Complete",
#                    stringsAsFactors = FALSE),
#         res$summary_row)
# }) %>%
#   dplyr::bind_rows()

# wb_sum <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb_sum, "PLS-DA_Summary")
# openxlsx::writeData(wb_sum, "PLS-DA_Summary", summary_all)
# openxlsx::saveWorkbook(wb_sum,
#                         file.path(out, "plsda_summary_all.xlsx"),
#                         overwrite = TRUE)

# print(summary_all)
# message("Script 12 complete → ", out)