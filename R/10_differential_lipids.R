# R/10_differential_lipids.R
# -------------------------------------------------------
# Differential lipid analysis: APO and PPWR_E
#
# Test:   Wilcoxon rank-sum or t-test (config differential.test)
# FDR:    toggle via config differential.apply_fdr
#   true  → BH correction, results → results/07_differential/fdr_BH_q{cutoff}/
#   false → raw p-values,  results → results/07_differential/raw_p{cutoff}/
#
# Volcano plots:
#   Figure 4A → results/figures/panels/fig4/fig4_A.pdf  (APO)
#   Figure 5A → results/figures/panels/fig5/fig5_A.pdf  (PPWR)
# -------------------------------------------------------
source("R/00_load_packages.R")
source("R/00_figure_theme.R")
cfg <- yaml::read_yaml("config/config.yml")

# ── Settings ──────────────────────────────────────────────
diff_cfg    <- cfg$differential
outcomes    <- diff_cfg$outcomes
test_method <- diff_cfg$test
apply_fdr   <- isTRUE(diff_cfg$apply_fdr)
fdr_method  <- diff_cfg$fdr_method
fdr_cutoff  <- diff_cfg$fdr_cutoff
raw_cutoff  <- diff_cfg$raw_p_cutoff
sig_cutoff  <- if (apply_fdr) fdr_cutoff else raw_cutoff
top_n       <- diff_cfg$top_n_labels

# ── Output subfolder based on method ─────────────────────
subfolder <- if (apply_fdr) {
  paste0("fdr_BH_q", gsub("\\.", "", sprintf("%.2f", fdr_cutoff)))
} else {
  paste0("raw_p", gsub("\\.", "", sprintf("%.2f", raw_cutoff)))
}
out <- file.path(cfg$output$s10, subfolder)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

message("Outcomes: ", paste(outcomes, collapse = ", "))
message("Test: ", test_method)
message("Mode: ", if (apply_fdr) paste0("FDR (BH, q<", fdr_cutoff, ")") else paste0("Raw p<", raw_cutoff))
message("Output: ", out)

tp_names  <- c("Identified_Baseline","Identified_TP36_38weeks","Identified_Postpartum")
tp_labels <- c("Baseline","Wk36_38","Postpartum")
tp_display <- c("10\u201316 weeks","36\u201338 weeks","Postpartum")

multiExpr_ID <- readRDS(file.path(cfg$output$processed, "multiExpr_ID.rds"))
demographic  <- readRDS(file.path(cfg$output$processed, "demographic_data_bmi.rds"))

# ── Test function ─────────────────────────────────────────
run_test <- function(expr, group_vec, outcome_name, tp_label) {
  common    <- intersect(rownames(expr), names(group_vec))
  expr      <- expr[common, ]
  group_vec <- group_vec[common]

  results <- lapply(colnames(expr), function(lipid) {
    vals  <- expr[[lipid]]
    grp   <- group_vec
    valid <- complete.cases(vals, grp)
    vals  <- vals[valid]; grp <- grp[valid]
    if (length(unique(grp)) < 2 || sum(valid) < 6) return(NULL)

    g0 <- vals[grp == 0]; g1 <- vals[grp == 1]

    tst <- if (test_method == "wilcoxon") {
      tryCatch(wilcox.test(g0, g1), error = function(e) NULL)
    } else {
      tryCatch(t.test(g0, g1), error = function(e) NULL)
    }
    if (is.null(tst)) return(NULL)

    data.frame(
      Lipid     = lipid,
      Outcome   = outcome_name,
      Timepoint = tp_label,
      log2FC    = mean(g1, na.rm = TRUE) - mean(g0, na.rm = TRUE),
      p_value   = tst$p.value,
      mean_ctrl = mean(g0, na.rm = TRUE),
      mean_case = mean(g1, na.rm = TRUE),
      n_ctrl    = sum(grp == 0),
      n_case    = sum(grp == 1),
      row.names = NULL
    )
  }) %>% bind_rows()

  results[order(results$p_value), ]
}

# ── Run tests ─────────────────────────────────────────────
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
    tp   <- tp_labels[i]
    expr <- multiExpr_ID[[tp_names[i]]]$data
    res  <- run_test(expr, grp_vec, outcome, tp)
    if (is.null(res) || nrow(res) == 0) next

    if (apply_fdr) {
      res$q_value_BH  <- p.adjust(res$p_value, method = fdr_method)
      res$significant <- res$q_value_BH < fdr_cutoff
    } else {
      res$q_value_BH  <- NA_real_
      res$significant <- res$p_value < raw_cutoff
    }
    res$direction <- ifelse(res$log2FC > 0, "UP", "DOWN")

    key <- paste0(tp, "_", outcome)
    diff_results_all[[key]] <- res
    write.xlsx(res, file.path(out, paste0("differential_", key, ".xlsx")))

    message(sprintf("  %s x %s: %d/%d significant (%s<%s)",
                    tp, outcome,
                    sum(res$significant, na.rm = TRUE), nrow(res),
                    if (apply_fdr) "q" else "p", sig_cutoff))
  }
}

# Save RDS with subfolder tag so 08 knows which version was used
saveRDS(diff_results_all,
        file.path(cfg$output$processed,
                  paste0("differential_lipids_", subfolder, ".rds")))
# Also save as generic for downstream compatibility
saveRDS(diff_results_all,
        file.path(cfg$output$processed, "differential_lipids_all.rds"))

# ── Volcano plots ─────────────────────────────────────────
make_volcano <- function(outcome, fig_id, panel_id, pos_label, neg_label) {
  keys <- paste0(tp_labels, "_", outcome)
  keys <- keys[keys %in% names(diff_results_all)]
  if (length(keys) == 0) {
    message("  No results for '", outcome, "' — skipping volcano"); return(invisible(NULL))
  }

  plots <- lapply(seq_along(keys), function(i) {
    res  <- diff_results_all[[keys[i]]]
    tp_d <- tp_display[i]
    if (is.null(res) || nrow(res) == 0) return(ggplot() + theme_void())

    top_sig <- res %>%
      filter(significant) %>%
      arrange(p_value) %>%
      slice_head(n = top_n)

    res$label <- ifelse(res$Lipid %in% top_sig$Lipid, res$Lipid, "")
    res$color  <- case_when(
      res$significant & res$log2FC > 0 ~ "UP",
      res$significant & res$log2FC < 0 ~ "DOWN",
      TRUE ~ "NS"
    )

    sig_line <- if (apply_fdr && any(res$significant, na.rm = TRUE)) {
      -log10(max(res$p_value[res$significant], na.rm = TRUE))
    } else if (!apply_fdr) {
      -log10(raw_cutoff)
    } else NA

    p <- ggplot(res, aes(x = log2FC, y = -log10(p_value),
                         colour = color, label = label)) +
      geom_point(size = 0.9, alpha = 0.7) +
      scale_colour_manual(
        values = c("UP" = "#E63946", "DOWN" = "#1D3D8F", "NS" = "#CCCCCC"),
        labels = c("UP" = pos_label, "DOWN" = neg_label, "NS" = "ns"),
        name   = NULL
      ) +
      labs(
        title    = tp_d,
        subtitle = paste0(
          sum(res$significant, na.rm = TRUE), " significant (",
          if (apply_fdr) paste0("q<", fdr_cutoff) else paste0("p<", raw_cutoff), ")"
        ),
        x = "log\u2082 fold change",
        y = "-log\u2081\u2080(p)"
      ) +
      theme_growell(grid = "both") +
      theme(legend.position = "bottom",
            legend.text     = element_text(size = FS_SMALL),
            plot.title      = element_text(size = FS_BASE, face = "bold"),
            plot.subtitle   = element_text(size = FS_SMALL,
                                           colour = GROWELL_COLORS$grey_dark))

    if (!is.na(sig_line))
      p <- p + geom_hline(yintercept = sig_line, linetype = "dashed",
                           colour = GROWELL_COLORS$grey_dark, linewidth = 0.3)

    if (requireNamespace("ggrepel", quietly = TRUE) && nrow(top_sig) > 0)
      p <- p + ggrepel::geom_text_repel(size = FS_SMALL / ggplot2::.pt,
                                         max.overlaps = 20, segment.size = 0.2,
                                         show.legend = FALSE)
    p
  })

  combined <- wrap_plots(plots, nrow = 1) +
    plot_annotation(
      theme = theme(plot.title = element_text(size = FS_TITLE, face = "bold",
                                               family = FONT_FAMILY))
    )

  ggsave(file.path(out, paste0("volcano_", outcome, ".pdf")),
         plot = combined, width = FIG_WIDTH_FULL, height = 80, units = "mm",
         device = cairo_pdf)
  ggsave(file.path(out, paste0("volcano_", outcome, ".png")),
         plot = combined, width = FIG_WIDTH_FULL, height = 80, units = "mm",
         dpi = FIGURE_DPI)

  save_panel(combined, fig_id, panel_id,
             width_mm = FIG_WIDTH_FULL, height_mm = 80)

  message("  Volcano: ", fig_id, "_", panel_id, " (", outcome, ")")
  invisible(combined)
}

message("\nGenerating volcano plots...")
make_volcano("apo",    "fig4", "A", "Higher in APO",  "Lower in APO")
make_volcano("ppwr_e", "fig5", "A", "Higher in PPWR", "Lower in PPWR")

# ── Summary table ─────────────────────────────────────────
summary_rows <- lapply(names(diff_results_all), function(key) {
  res <- diff_results_all[[key]]
  if (is.null(res) || nrow(res) == 0) return(NULL)
  data.frame(
    Analysis      = key,
    Total         = nrow(res),
    Sig_raw_p05   = sum(res$p_value < 0.05, na.rm = TRUE),
    Sig_threshold = sum(res$significant, na.rm = TRUE),
    N_UP          = sum(res$significant & res$direction == "UP",   na.rm = TRUE),
    N_DOWN        = sum(res$significant & res$direction == "DOWN", na.rm = TRUE)
  )
}) %>% bind_rows()

write.xlsx(summary_rows, file.path(out, "differential_summary.xlsx"))
print(summary_rows)

message("\nScript 10 complete → ", out)
message("Mode: ", subfolder)
message("To switch: change differential.apply_fdr in config/config.yml")