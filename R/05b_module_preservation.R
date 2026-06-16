# R/05b_module_preservation.R
# -------------------------------------------------------
# Module preservation analysis with Z summary statistics
#
# Z summary interpretation (Langfelder & Horvath 2011):
#   Z summary > 10  → strongly preserved
#   Z summary 2-10  → moderately preserved
#   Z summary < 2   → not preserved
#
# Data structure confirmed from diagnostic:
#   Zsummary    → mp$preservation$Z[[1]][[2]]
#   medianRank  → mp$preservation$observed[[1]][[2]]
#   moduleSize  → in BOTH dataframes as first column
#
# Outputs → results/05_consensus_wgcna/<run_tag>/
# -------------------------------------------------------

source("R/00_load_packages.R")
source("R/00_figure_theme.R")

library(ggplot2)
library(dplyr)

`%||%` <- function(a, b) if (!is.null(a)) a else b

cfg <- yaml::read_yaml("config/config.yml")

# ── Read run tag ──────────────────────────────────────────
wgcna_run_tag <- trimws(readLines(
  file.path(cfg$output$processed, "current_wgcna_run_tag.txt")))
message("WGCNA run: ", wgcna_run_tag)

out <- file.path(cfg$output$s05, wgcna_run_tag)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ── Load data ─────────────────────────────────────────────
multiExpr_ID     <- readRDS(file.path(cfg$output$processed,
                                       "multiExpr_ID.rds"))
consensusMods_ID <- readRDS(file.path(cfg$output$processed,
  paste0("consensusMods_ID_", wgcna_run_tag, ".rds")))

module_colors <- consensusMods_ID$colors
message("Modules to test: ",
        length(unique(module_colors[module_colors != "grey"])),
        " non-grey modules")

# ── Timepoint definitions ─────────────────────────────────
tp_names   <- c("Identified_Baseline",
                "Identified_TP36_38weeks",
                "Identified_Postpartum")
tp_short   <- c("Baseline", "Wk36_38", "Postpartum")
tp_display <- c("10-16 weeks", "36-38 weeks", "Postpartum")

# ── Pairwise comparisons ──────────────────────────────────
comparisons <- list(
  list(ref = 1, test = 2,
       label   = paste0(tp_short[1],   "_vs_", tp_short[2]),
       display = paste0(tp_display[1], " vs ", tp_display[2])),
  list(ref = 1, test = 3,
       label   = paste0(tp_short[1],   "_vs_", tp_short[3]),
       display = paste0(tp_display[1], " vs ", tp_display[3])),
  list(ref = 2, test = 3,
       label   = paste0(tp_short[2],   "_vs_", tp_short[3]),
       display = paste0(tp_display[2], " vs ", tp_display[3]))
)

# ── Helper: extract stats from mp object ──────────────────
# Confirmed structure from diagnostic:
#   Z stats    → mp$preservation$Z[[1]][[2]]
#   Obs stats  → mp$preservation$observed[[1]][[2]]
extract_mp_stats <- function(mp) {

  z_df  <- mp$preservation$Z[[1]][[2]]
  obs_df <- mp$preservation$observed[[1]][[2]]

  # Both dataframes are indexed by module name (rownames)
  # Merge on rownames so columns align correctly
  shared_modules <- intersect(rownames(z_df), rownames(obs_df))

  combined <- data.frame(
    Module       = shared_modules,
    moduleSize   = z_df[shared_modules,   "moduleSize"],
    Zsummary     = round(z_df[shared_modules,   "Zsummary.pres"],    2),
    Zdensity     = round(z_df[shared_modules,   "Zdensity.pres"],    2),
    Zconnectivity= round(z_df[shared_modules,   "Zconnectivity.pres"], 2),
    medianRank   = round(obs_df[shared_modules, "medianRank.pres"],  2),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  combined
}

# ── Run preservation for each comparison ──────────────────
message("\nRunning module preservation analyses...")
message("nPermutations = 200 | quickCor = 1")
message("Each comparison takes ~15-30 min on a laptop\n")

preservation_results <- list()

for (comp in comparisons) {

  rds_path <- file.path(out, paste0("mp_", comp$label, ".rds"))

  # Load cached result if available — avoids rerunning permutations
  if (file.exists(rds_path)) {
    message(sprintf("  Loading cached: %s", comp$label))
    mp <- readRDS(rds_path)
  } else {
    message(sprintf("  Running: %s vs %s...", tp_short[comp$ref],
                    tp_short[comp$test]))

    multiData_pair <- list(
      list(data = multiExpr_ID[[tp_names[comp$ref]]]$data),
      list(data = multiExpr_ID[[tp_names[comp$test]]]$data)
    )

    mp <- WGCNA::modulePreservation(
      multiData         = multiData_pair,
      multiColor        = list(module_colors, module_colors),
      referenceNetworks = 1,
      nPermutations     = 200,
      randomSeed        = 12345,
      quickCor          = 1,
      verbose           = 3
    )

    # Save immediately after completion — prevents losing work if
    # later sections crash
    saveRDS(mp, rds_path)
    message(sprintf("  Saved: mp_%s.rds", comp$label))
  }

  # Extract stats using confirmed column locations
  stats_df <- extract_mp_stats(mp)

  preservation_results[[comp$label]] <- list(
    mp       = mp,
    stats    = stats_df,
    label    = comp$label,
    display  = comp$display
  )

  # ── Standard WGCNA preservation plots ──────────────────
  stats_plot <- stats_df[stats_df$Module != "grey", ]

  pdf(file.path(out,
    paste0("module_preservation_", comp$label, ".pdf")),
    width = 10, height = 5)
  par(mfrow = c(1, 2))
  par(mar = c(4.5, 4.5, 2.5, 1))

  # Z summary vs module size
  plot(stats_plot$moduleSize, stats_plot$Zsummary,
       col  = stats_plot$Module,
       pch  = 19, cex = 1.4,
       xlab = "Module size (number of lipids)",
       ylab = "Z summary",
       main = paste0("Preservation: ", comp$display),
       ylim = c(min(-2, min(stats_plot$Zsummary, na.rm = TRUE)),
                max(15, max(stats_plot$Zsummary, na.rm = TRUE))))
  WGCNA::labelPoints(
    stats_plot$moduleSize, stats_plot$Zsummary,
    stats_plot$Module, cex = 0.6, offs = 0.04)
  abline(h = 2,  col = "blue",  lty = 2, lwd = 1.5)
  abline(h = 10, col = "green", lty = 2, lwd = 1.5)
  legend("topright",
         legend = c("Z > 10: strongly preserved",
                    "Z 2-10: moderately preserved",
                    "Z < 2: not preserved"),
         col = c("green", "blue", "red"),
         lty = 2, lwd = 1.5, cex = 0.7, bg = "white")

  # Median rank vs module size
  plot(stats_plot$moduleSize, stats_plot$medianRank,
       col  = stats_plot$Module,
       pch  = 19, cex = 1.4,
       xlab = "Module size (number of lipids)",
       ylab = "Median rank",
       main = paste0("Median rank: ", comp$display))
  WGCNA::labelPoints(
    stats_plot$moduleSize, stats_plot$medianRank,
    stats_plot$Module, cex = 0.6, offs = 0.04)

  dev.off()
  message(sprintf("  Plot saved: module_preservation_%s.pdf", comp$label))
}

# ── Aggregate all comparisons ─────────────────────────────
message("\nAggregating Z summary statistics...")

zsummary_all <- lapply(names(preservation_results), function(nm) {
  res <- preservation_results[[nm]]
  df  <- res$stats
  df$Comparison <- res$display
  df$Preservation <- dplyr::case_when(
  df$Module %in% c("grey", "gold") ~ "Internal control",
  df$Zsummary >= 10  ~ "Strong (Z >= 10)",
  df$Zsummary >= 2   ~ "Moderate (2 <= Z < 10)",
  TRUE               ~ "Not preserved (Z < 2)"
  )
  df
}) %>% dplyr::bind_rows()

# ── Save Excel ────────────────────────────────────────────
wb <- openxlsx::createWorkbook()

# Sheet 1: full long-format table
openxlsx::addWorksheet(wb, "Z_Summary_Long")
openxlsx::writeData(wb, "Z_Summary_Long", zsummary_all)

# Colour coding
strong_style   <- openxlsx::createStyle(
  fgFill = "#C8E6C9", fontColour = "#1B5E20", textDecoration = "bold")
moderate_style <- openxlsx::createStyle(
  fgFill = "#FFF9C4", fontColour = "#F57F17")
weak_style     <- openxlsx::createStyle(
  fgFill = "#FFCDD2", fontColour = "#B71C1C")

for (i in seq_len(nrow(zsummary_all))) {
  style <- if (grepl("Strong",   zsummary_all$Preservation[i])) strong_style   else
           if (grepl("Moderate", zsummary_all$Preservation[i])) moderate_style else
           weak_style
  openxlsx::addStyle(wb, "Z_Summary_Long", style,
                      rows = i + 1, cols = 1:ncol(zsummary_all))
}

# Sheet 2: wide format — one row per module
zsummary_wide <- zsummary_all %>%
   filter(!Module %in% c("grey", "gold")) %>%
  dplyr::select(Module, moduleSize, Comparison, Zsummary, Preservation) %>%
  tidyr::pivot_wider(
    names_from  = Comparison,
    values_from = c(Zsummary, Preservation),
    names_sep   = " | "
  )

openxlsx::addWorksheet(wb, "Z_Summary_Wide")
openxlsx::writeData(wb, "Z_Summary_Wide", zsummary_wide)
openxlsx::setColWidths(wb, "Z_Summary_Wide",
                        cols = 1:ncol(zsummary_wide), widths = 22)

openxlsx::saveWorkbook(wb,
  file.path(out, "module_preservation_summary.xlsx"),
  overwrite = TRUE)
message("  Saved: module_preservation_summary.xlsx")

# ── Reset graphics device state before ggplot2 figures ───
# Base R pdf() calls earlier in the script may leave device state
# that conflicts with ggsave + ggrepel rendering
while (!is.null(dev.list())) dev.off()

# ── Figure S6: publication quality ────────────────────────
# Replace the Figure S6 save section in 05b_module_preservation.R
# with this approach — saves panels separately, avoids patchwork viewport error

# ── Figure S6: publication quality ────────────────────────
message("\nGenerating Figure S6...")

plot_df <- zsummary_all %>%
  filter(!Module %in% c("grey", "gold")) %>%
  mutate(
    Preservation = factor(Preservation, levels = c(
      "Strong (Z >= 10)",
      "Moderate (2 <= Z < 10)",
      "Not preserved (Z < 2)"))
  )

# Module dot colours
safe_color <- function(col_name) {
  tryCatch(
    { grDevices::col2rgb(col_name); col_name },
    error = function(e) "grey50"
  )
}
dot_colors <- setNames(
  sapply(unique(plot_df$Module), safe_color),
  unique(plot_df$Module)
)
dot_colors["white"] <- "grey90"

# Shared theme
theme_preservation <- function() {
  ggplot2::theme_minimal(base_family = "Helvetica", base_size = 8) +
  ggplot2::theme(
    panel.border     = ggplot2::element_rect(
                         colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = ggplot2::element_line(
                         colour = "#EEEEEE", linewidth = 0.2),
    panel.grid.minor = ggplot2::element_blank(),
    strip.text       = ggplot2::element_text(
                         size = 7, face = "bold", family = "Helvetica"),
    axis.text        = ggplot2::element_text(
                         size = 7, face = "bold",
                         colour = "black", family = "Helvetica"),
    axis.title       = ggplot2::element_text(
                         size = 8, face = "bold", family = "Helvetica"),
    plot.title       = ggplot2::element_text(
                         size = 9, face = "bold", family = "Helvetica"),
    plot.tag         = ggplot2::element_text(
                         size = 10, face = "bold", family = "Helvetica"),
    plot.margin      = ggplot2::margin(4, 6, 4, 4, unit = "mm"),
    legend.position  = "none"
  )
}

# Panel a: Z summary vs module size
p_zsummary <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(x = moduleSize, y = Zsummary,
               colour = Module, label = Module)
) +
  ggplot2::geom_hline(yintercept = 2,  linetype = "dashed",
                       colour = "#1565C0", linewidth = 0.5) +
  ggplot2::geom_hline(yintercept = 10, linetype = "dashed",
                       colour = "#2E7D32", linewidth = 0.5) +
  ggplot2::geom_point(size = 2.5, alpha = 0.9) +
  ggplot2::geom_text(nudge_x      = 0.5,
                    nudge_y      = 0.3,
                    size         = 2.0,
                    family       = "Helvetica",
                    check_overlap = TRUE,
                    show.legend  = FALSE
                    ) +
  ggplot2::annotate("text", x = Inf, y = 10.4,
                     label = "Z = 10 (strongly preserved)",
                     hjust = 1.05, size = 2.2,
                     colour = "#2E7D32", family = "Helvetica") +
  ggplot2::annotate("text", x = Inf, y = 2.4,
                     label = "Z = 2 (moderately preserved)",
                     hjust = 1.05, size = 2.2,
                     colour = "#1565C0", family = "Helvetica") +
  ggplot2::scale_colour_manual(values = dot_colors, guide = "none") +
  ggplot2::facet_wrap(~ Comparison, nrow = 1) +
  ggplot2::labs(
    tag   = "a",
    x     = "Module size (number of lipids)",
    y     = "Z summary",
    title = "Module preservation across timepoints"
  ) +
  theme_preservation()

# Panel b: Median rank vs module size
p_medrank <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(x = moduleSize, y = medianRank,
               colour = Module, label = Module)
) +
  ggplot2::geom_point(size = 2.5, alpha = 0.9) +
  ggplot2::geom_text(nudge_x      = 0.5,
                    nudge_y      = 0.3,
                    size         = 2.0,
                    family       = "Helvetica",
                    check_overlap = TRUE,
                    show.legend  = FALSE
                    ) +
  ggplot2::scale_colour_manual(values = dot_colors, guide = "none") +
  ggplot2::facet_wrap(~ Comparison, nrow = 1) +
  ggplot2::labs(
    tag   = "b",
    x     = "Module size (number of lipids)",
    y     = "Median rank",
    title = "Median rank (lower = better preserved)"
  ) +
  theme_preservation()

# ── Save panels separately — avoids patchwork viewport error ──
fig_s6_base <- file.path(out, "FigureS6_module_preservation")

# Save panel a
ggplot2::ggsave(
  filename = paste0(fig_s6_base, "_a_Zsummary.pdf"),
  plot     = p_zsummary,
  width    = FIG_WIDTH_FULL,
  height   = 85,
  units    = "mm",
  dpi      = 600,
  device   = grDevices::cairo_pdf
)

ggplot2::ggsave(
  filename = paste0(fig_s6_base, "_a_Zsummary.tiff"),
  plot     = p_zsummary,
  width    = FIG_WIDTH_FULL,
  height   = 85,
  units    = "mm",
  dpi      = 600,
  device   = function(filename, ...) grDevices::tiff(
    filename, ..., compression = "lzw",
    type = "cairo", family = "Helvetica")
)

# Save panel b
ggplot2::ggsave(
  filename = paste0(fig_s6_base, "_b_MedianRank.pdf"),
  plot     = p_medrank,
  width    = FIG_WIDTH_FULL,
  height   = 85,
  units    = "mm",
  dpi      = 600,
  device   = grDevices::cairo_pdf
)

ggplot2::ggsave(
  filename = paste0(fig_s6_base, "_b_MedianRank.tiff"),
  plot     = p_medrank,
  width    = FIG_WIDTH_FULL,
  height   = 85,
  units    = "mm",
  dpi      = 600,
  device   = function(filename, ...) grDevices::tiff(
    filename, ..., compression = "lzw",
    type = "cairo", family = "Helvetica")
)

message("Figure S6 panels saved separately:")
message("  a: ", fig_s6_base, "_a_Zsummary.pdf/.tiff")
message("  b: ", fig_s6_base, "_b_MedianRank.pdf/.tiff")
message("  Combine panels a and b in Illustrator for final Figure S6")



# ── Console summary ───────────────────────────────────────
message("\n── Z summary statistics summary ──")
for (comp_nm in unique(zsummary_all$Comparison)) {
  sub <- zsummary_all %>%
    filter(Comparison == comp_nm, Module != "grey")
  message(sprintf("\n  %s:", comp_nm))
  message(sprintf("    Strongly preserved  (Z >= 10): %d modules",
                  sum(grepl("Strong",   sub$Preservation))))
  message(sprintf("    Moderately preserved (2-10):   %d modules",
                  sum(grepl("Moderate", sub$Preservation))))
  message(sprintf("    Not preserved (Z < 2):         %d modules",
                  sum(grepl("Not",      sub$Preservation))))
}

message("\nRebuttal note (WM-4):")
message("  Figure S6 now includes Z summary and median rank statistics.")
message("  Z >= 10: strongly preserved | Z 2-10: moderately preserved")
message("\nScript 05b complete -> ", out)