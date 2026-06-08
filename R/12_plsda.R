# R/12_plsda.R
# -------------------------------------------------------
# PLS-DA with ropls backend
#   - Q², R²X / R²Y, permutation p-value (native ropls)
#   - Biplot (scores + loadings overlay, ggplot2 + ggrepel)
#   - Config-driven outcome selection via plsda.outcomes
#   - Per-timepoint × outcome Excel summary
# Outputs → results/09_plsda/
# -------------------------------------------------------

source("R/00_load_packages.R")
source("R/00_figure_theme.R")

cfg        <- yaml::read_yaml("config/config.yml")
out        <- cfg$output$s12
dir.create(out, showWarnings = FALSE, recursive = TRUE)

multiExpr_ID  <- readRDS(file.path(cfg$output$processed, "multiExpr_ID.rds"))
status_labels <- readRDS(file.path(cfg$output$processed, "lipid_status_labels.rds"))

tp_names  <- names(multiExpr_ID)
tp_labels <- c("Baseline", "Wk36_38", "Postpartum")

# ── Outcomes: pull from config, fall back to all available ──────────────────
cfg_outcomes <- cfg$plsda$outcomes          # character vector or NULL
all_outcomes <- colnames(status_labels)
outcomes <- if (!is.null(cfg_outcomes) && length(cfg_outcomes) > 0) {
  missing_oc <- setdiff(cfg_outcomes, all_outcomes)
  if (length(missing_oc) > 0)
    warning("plsda.outcomes not found in status_labels: ",
            paste(missing_oc, collapse = ", "))
  intersect(cfg_outcomes, all_outcomes)
} else {
  all_outcomes
}
message("Running PLS-DA for outcomes: ", paste(outcomes, collapse = ", "))

# ── PLS-DA parameters from config ───────────────────────────────────────────
ncomp   <- cfg$plsda$ncomp      # number of components
permI   <- cfg$plsda$permI      # permutation iterations (e.g. 999)
crossI  <- cfg$plsda$crossI     # cross-validation folds (e.g. 7)
n_top   <- cfg$plsda$n_top_loadings %||% 10L   # top loadings to show in biplot

# ── Colour palette (group labels come from factor levels) ───────────────────
group_colours <- cfg$plsda$group_colours   # named list from config, or NULL
default_pal   <- c(
  "#4DAFCA", "#E94D36", "#4DAF4A", "#984EA3",
  "#FF7F00", "#A65628", "#F781BF", "#999999"
)

# ── Helper: %||% (null-coalescing) ───────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Biplot function ──────────────────────────────────────────────────────────
#' @param fit        ropls opls object (2-component PLS-DA)
#' @param group_vec  named factor/character; names = sample IDs
#' @param tp_label   string for plot title
#' @param outcome_label string for plot subtitle
#' @param n_top      number of top-VIP loadings to display
#' @param pal        named character vector of colours keyed to group levels,
#'                   or NULL to use defaults
make_plsda_biplot <- function(fit, group_vec, tp_label, outcome_label,
                               n_top = 10, pal = NULL) {

  # --- scores -----------------------------------------------------------------
  sc <- as.data.frame(ropls::getScoreMN(fit))
  if (ncol(sc) < 2) {
    message("  Biplot skipped: fewer than 2 components in model")
    return(NULL)
  }
  colnames(sc)[1:2] <- c("Comp1", "Comp2")
  sc$Sample <- rownames(sc)
  sc$Group  <- as.character(group_vec[rownames(sc)])

  # --- loadings (X-side) ------------------------------------------------------
  ld <- as.data.frame(ropls::getLoadingMN(fit))
  colnames(ld)[1:2] <- c("Comp1", "Comp2")
  ld$Lipid <- rownames(ld)

  # top n by VIP on Comp1
  vip_vals  <- ropls::getVipVn(fit)
  top_names <- names(sort(vip_vals, decreasing = TRUE))[seq_len(min(n_top, nrow(ld)))]
  ld_top    <- ld[ld$Lipid %in% top_names, ]

  # --- axis labels (R²X per component, matching MetaboAnalyst convention) ----
  summ <- ropls::getSummaryDF(fit)
  # modelDF rows correspond to each component; R2X(cum) is cumulative
  # Back-calculate per-component R²X
  r2x_cum <- fit@modelDF[["R2X(cum)"]]
  r2x_per <- diff(c(0, r2x_cum))
  xlab <- sprintf("Component 1 (%.1f%%)", r2x_per[1] * 100)
  ylab <- if (length(r2x_per) >= 2)
    sprintf("Component 2 (%.1f%%)", r2x_per[2] * 100)
  else
    "Component 2"

  # --- scale arrows into score space ------------------------------------------
  score_range <- max(abs(c(range(sc$Comp1), range(sc$Comp2))), na.rm = TRUE)
  load_range  <- max(abs(c(range(ld_top$Comp1), range(ld_top$Comp2))), na.rm = TRUE)
  arrow_scale <- if (load_range > 0) 0.75 * score_range / load_range else 1
  ld_top <- ld_top %>%
    dplyr::mutate(x_end = Comp1 * arrow_scale,
                  y_end = Comp2 * arrow_scale)

  # --- colour palette ---------------------------------------------------------
  grp_levels <- sort(unique(sc$Group))
  if (!is.null(pal) && all(grp_levels %in% names(pal))) {
    col_vals <- pal[grp_levels]
  } else {
    col_vals <- setNames(default_pal[seq_along(grp_levels)], grp_levels)
  }

  # --- plot -------------------------------------------------------------------
  ggplot2::ggplot() +
    # sample scores
    ggplot2::geom_point(
      data = sc,
      ggplot2::aes(x = Comp1, y = Comp2, colour = Group),
      size = 2.5, alpha = 0.85
    ) +
    # loading arrows
    ggplot2::geom_segment(
      data = ld_top,
      ggplot2::aes(x = 0, y = 0, xend = x_end, yend = y_end),
      arrow     = ggplot2::arrow(length = ggplot2::unit(0.18, "cm"),
                                  type = "closed"),
      colour    = "grey30",
      linewidth = 0.4
    ) +
    # lipid labels — ggrepel avoids overlap
    ggrepel::geom_text_repel(
      data = ld_top,
      ggplot2::aes(x = x_end, y = y_end, label = Lipid),
      size         = 2.8,
      colour       = "grey15",
      max.overlaps = 25,
      segment.size = 0.2,
      box.padding  = 0.35
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                         colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                         colour = "grey70", linewidth = 0.3) +
    ggplot2::scale_colour_manual(values = col_vals, name = NULL) +
    ggplot2::labs(
      x        = xlab,
      y        = ylab,
      title    = paste("PLS-DA Biplot:", tp_label),
      subtitle = outcome_label
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position  = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle    = ggplot2::element_text(size = 9, colour = "grey40")
    )
}

# ── Core fitting function ────────────────────────────────────────────────────
#' Fit a PLS-DA model with ropls, export plots + Excel, return summary row
#'
#' @return list with slots: model, summary_row, vip_df, or NULL on failure
run_plsda_ropls <- function(expr, group_vec, outcome_name, tp_label) {

  # --- align samples ----------------------------------------------------------
  common <- intersect(rownames(expr), names(group_vec))
  if (length(common) == 0) {
    message("  Skipping ", tp_label, " / ", outcome_name, ": no overlapping samples")
    return(NULL)
  }
  X <- expr[common, , drop = FALSE]
  Y <- factor(group_vec[common])

  # --- filter constant / all-NA columns ---------------------------------------
  valid_cols <- apply(X, 2, function(x) !all(is.na(x)) && var(x, na.rm = TRUE) > 0)
  X <- X[, valid_cols, drop = FALSE]

  # --- complete cases ---------------------------------------------------------
  cc <- stats::complete.cases(X, as.integer(Y))
  X <- X[cc, , drop = FALSE]
  Y <- Y[cc]

  if (length(unique(Y)) < 2) {
    message("  Skipping ", tp_label, " / ", outcome_name,
            ": outcome has < 2 levels after filtering")
    return(NULL)
  }
  if (nrow(X) < 10) {
    message("  Skipping ", tp_label, " / ", outcome_name,
            ": n = ", nrow(X), " < 10 after filtering")
    return(NULL)
  }

  message(sprintf("  PLS-DA [ropls]: %s × %s  (n=%d, p=%d, permI=%d)",
                  tp_label, outcome_name, nrow(X), ncol(X), permI))

  # --- fit ropls opls ---------------------------------------------------------
  fit <- tryCatch(
  ropls::opls(
    x     = X,
    y     = Y,
    predI = ncomp,
    permI = permI,
    crossvalI = crossI,        # was crossI in older versions
    info.txtC = "none",        # suppresses printed output
    fig.pdfC  = "none"         # suppresses auto PDF generation
  ),
  error = function(e) {
    message("  ropls error: ", conditionMessage(e))
    NULL
  }
)
  if (is.null(fit)) return(NULL)

  # --- extract summary stats --------------------------------------------------
  mdf      <- fit@modelDF                     # one row per component
  r2x_cum  <- mdf[["R2X(cum)"]]
  r2y_cum  <- mdf[["R2Y(cum)"]]
  q2_cum   <- mdf[["Q2(cum)"]]
  r2x_per  <- diff(c(0, r2x_cum))
  r2y_per  <- diff(c(0, r2y_cum))

  # permutation p-values live in summaryDF (one row per perm stat)
  sdf      <- fit@summaryDF
  perm_pR2Y <- if ("pR2Y" %in% colnames(sdf)) sdf[nrow(sdf), "pR2Y"] else NA_real_
  perm_pQ2  <- if ("pQ2"  %in% colnames(sdf)) sdf[nrow(sdf), "pQ2"]  else NA_real_

  # VIP scores (single vector, Comp1-focused from ropls)
  vip_vals <- ropls::getVipVn(fit)
  vip_df   <- data.frame(
    Lipid     = names(vip_vals),
    VIP_Comp1 = vip_vals,
    row.names = NULL
  ) %>% dplyr::arrange(dplyr::desc(VIP_Comp1))

  # --- output prefix ----------------------------------------------------------
  pfx <- file.path(out, paste0("plsda_", tp_label, "_", outcome_name))

  # --- biplot (publication quality TIFF + PDF) --------------------------------
  bp <- make_plsda_biplot(
    fit            = fit,
    group_vec      = group_vec,
    tp_label       = tp_label,
    outcome_label  = outcome_name,
    n_top          = n_top,
    pal            = group_colours
  )
  if (!is.null(bp)) {
    ggplot2::ggsave(paste0(pfx, "_biplot.pdf"),  bp, width = 7,   height = 6,   dpi = 300)
    ggplot2::ggsave(paste0(pfx, "_biplot.tiff"), bp, width = 7,   height = 6,   dpi = 600,
                    compression = "lzw")
  }

  # --- permutation plot (ropls built-in, saved to PDF) ------------------------
  pdf(paste0(pfx, "_permutation.pdf"), width = 7, height = 5)
  tryCatch(ropls::plot(fit, typeVc = "permutation"), error = function(e) NULL)
  dev.off()

  # --- score plot (ropls built-in) --------------------------------------------
  pdf(paste0(pfx, "_scores.pdf"), width = 7, height = 6)
  tryCatch(ropls::plot(fit, typeVc = "x-score"), error = function(e) NULL)
  dev.off()

  # --- Excel output -----------------------------------------------------------
  wb <- openxlsx::createWorkbook()

  # Sheet 1: model summary
  n_comp_fit <- nrow(mdf)
  summary_sheet <- data.frame(
    Timepoint   = tp_label,
    Outcome     = outcome_name,
    N_samples   = nrow(X),
    N_lipids    = ncol(X),
    N_comp      = n_comp_fit,
    CrossI      = crossI,
    PermI       = permI
  )
  for (k in seq_len(n_comp_fit)) {
    summary_sheet[[paste0("R2X_Comp",  k)]] <- round(r2x_per[k],  4)
    summary_sheet[[paste0("R2Y_Comp",  k)]] <- round(r2y_per[k],  4)
    summary_sheet[[paste0("Q2_cum_Comp", k)]] <- round(q2_cum[k], 4)
  }
  summary_sheet$perm_pR2Y <- perm_pR2Y
  summary_sheet$perm_pQ2  <- perm_pQ2

  openxlsx::addWorksheet(wb, "Model_Summary")
  openxlsx::writeData(wb, "Model_Summary", summary_sheet)

  # Sheet 2: full modelDF from ropls
  openxlsx::addWorksheet(wb, "modelDF")
  openxlsx::writeData(wb, "modelDF", as.data.frame(mdf))

  # Sheet 3: VIP scores
  openxlsx::addWorksheet(wb, "VIP_scores")
  openxlsx::writeData(wb, "VIP_scores", vip_df)

  # Sheet 4: loadings
  ld_all <- as.data.frame(ropls::getLoadingMN(fit))
  colnames(ld_all) <- paste0("Comp", seq_len(ncol(ld_all)))
  ld_all$Lipid <- rownames(ld_all)
  ld_all <- ld_all[, c("Lipid", setdiff(names(ld_all), "Lipid"))]
  openxlsx::addWorksheet(wb, "Loadings")
  openxlsx::writeData(wb, "Loadings", ld_all)

  openxlsx::saveWorkbook(wb, paste0(pfx, ".xlsx"), overwrite = TRUE)

  # --- return -----------------------------------------------------------------
  list(
    model       = fit,
    vip_df      = vip_df,
    summary_row = summary_sheet,
    tp          = tp_label,
    outcome     = outcome_name,
    n_samples   = nrow(X),
    n_lipids    = ncol(X)
  )
}

# ── Main loop ────────────────────────────────────────────────────────────────
plsda_results <- list()

for (i in seq_along(tp_names)) {
  expr <- multiExpr_ID[[tp_names[i]]]$data
  tp   <- tp_labels[i]

  for (outcome in outcomes) {
    key     <- paste0(tp, "_", outcome)
    grp_vec <- setNames(status_labels[[outcome]], rownames(status_labels))

    plsda_results[[key]] <- run_plsda_ropls(
      expr         = expr,
      group_vec    = grp_vec,
      outcome_name = outcome,
      tp_label     = tp
    )
  }
}

saveRDS(plsda_results, file.path(cfg$output$processed, "plsda_results.rds"))

# ── Aggregate summary table ──────────────────────────────────────────────────
summary_all <- lapply(names(plsda_results), function(key) {
  res <- plsda_results[[key]]
  if (is.null(res)) {
    return(data.frame(Key = key, Status = "Failed or Skipped",
                      stringsAsFactors = FALSE))
  }
  cbind(data.frame(Key = key, Status = "Complete",
                   stringsAsFactors = FALSE),
        res$summary_row)
}) %>%
  dplyr::bind_rows()

wb_sum <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_sum, "PLS-DA_Summary")
openxlsx::writeData(wb_sum, "PLS-DA_Summary", summary_all)
openxlsx::saveWorkbook(wb_sum,
                        file.path(out, "plsda_summary_all.xlsx"),
                        overwrite = TRUE)

print(summary_all)
message("Script 12 complete → ", out)