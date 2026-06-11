# R/11_lipid_boxplots.R
# -------------------------------------------------------
# Boxplots of individual lipids of interest across 3 timepoints
# for APO and PPWR_E outcomes
#
# Figure 4C, 4D (APO) and Figure 5C, 5D (PPWR)
#
# Layout controlled by: plot_layout (in this script)
#   "vertical"   → 3 timepoints stacked top-to-bottom (ncol=1)
#                  good for side-by-side panels in a figure row
#   "horizontal" → 3 timepoints left-to-right (nrow=1)
#                  good for a single panel spanning full width
#
# Lipid names must match exactly as stored in multiExpr_ID
# (format: "TG 50:1", "TG 46:0" etc — run script to see examples)
# -------------------------------------------------------
source("R/00_load_packages.R")
source("R/00_figure_theme.R")
cfg <- yaml::read_yaml("config/config.yml")

out <- file.path(cfg$output$s11)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ── Layout toggle ─────────────────────────────────────────
# "vertical"   → timepoints stacked top-to-bottom  (89mm wide × 150mm tall)
# "horizontal" → timepoints left-to-right          (183mm wide × 60mm tall)
plot_layout <- "horizontal"   # ← change to "vertical" if needed

panel_width  <- if (plot_layout == "vertical") FIG_WIDTH_SINGLE else FIG_WIDTH_FULL
panel_height <- if (plot_layout == "vertical") 150 else 70
ncols        <- if (plot_layout == "vertical") 1 else 3

# ── Lipids of interest ────────────────────────────────────
# Use exact names from multiExpr_ID (spaces and colons)
# Check available names: grep('^TG', colnames(multiExpr_ID[[1]]$data), value=TRUE)
lipids_of_interest <- list(
  apo    = list(
    C = "TG 48:0",   # Figure 4C
    D = "TG 46:0"    # Figure 4D
  ),
  ppwr_e = list(
    C = "TG 58:9",   # Figure 5C
    D = "TG 58:8"    # Figure 5D
  )
)

# ── Outcome display info ──────────────────────────────────
outcome_info <- list(
  apo = list(
    no_label  = "No APO",
    yes_label = "APO",
    colors    = c("No APO"  = "#1D3D8F", "APO"  = "#E63946"),
    fig_id    = "fig4",
    xlab      = "APO status"
  ),
  ppwr_e = list(
    no_label  = "No PPWR",
    yes_label = "PPWR",
    colors    = c("No PPWR" = "#1D3D8F", "PPWR" = "#E63946"),
    fig_id    = "fig5",
    xlab      = "PPWR status"
  )
)

# ── Load data ─────────────────────────────────────────────
multiExpr_ID <- readRDS(file.path(cfg$output$processed, "multiExpr_ID.rds"))
demographic  <- readRDS(file.path(cfg$output$processed, "demographic_data_bmi.rds"))

tp_names   <- c("Identified_Baseline","Identified_TP36_38weeks","Identified_Postpartum")
tp_labels  <- c("Baseline","Wk36_38","Postpartum")
tp_display <- c("10\u201316 weeks","36\u201338 weeks","Postpartum")

message("Layout: ", plot_layout, " (", ncols, " column(s) per panel)")

# ── Single timepoint boxplot ──────────────────────────────
plot_one_tp <- function(lipid_name, tp_idx, outcome, info) {
  expr <- multiExpr_ID[[tp_names[tp_idx]]]$data
  tp_d <- tp_display[tp_idx]

  if (!(lipid_name %in% colnames(expr))) {
    return(ggplot() +
             annotate("text", x=0.5, y=0.5,
                      label=paste(tp_d, "\nnot available"),
                      size=FS_SMALL/ggplot2::.pt, colour="#999999") +
             theme_void() +
             labs(title=tp_d))
  }

  grp_vec <- setNames(demographic[[outcome]], rownames(demographic))
  common  <- intersect(rownames(expr), names(grp_vec))
  vals    <- expr[common, lipid_name]
  grps    <- grp_vec[common]

  df <- na.omit(data.frame(
    Intensity = as.numeric(vals),
    Group     = factor(grps, levels = c(0, 1),
                       labels = c(info$no_label, info$yes_label))
  ))

  tst   <- tryCatch(wilcox.test(Intensity ~ Group, data = df),
                    error = function(e) NULL)
  p_val <- if (!is.null(tst)) tst$p.value else NA
  stars <- ifelse(!is.na(p_val) & p_val < 0.001, "***",
           ifelse(!is.na(p_val) & p_val < 0.01,  "**",
           ifelse(!is.na(p_val) & p_val < 0.05,  "*", "ns")))
  subtitle <- if (!is.na(p_val)) paste0("p = ", signif(p_val, 2), " ", stars) else "p = NA"

  ggplot(df, aes(x = Group, y = Intensity, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7,
                 linewidth = 0.4, width = 0.55) +
    geom_jitter(aes(colour = Group), width = 0.12,
                size = 1, alpha = 0.6, show.legend = FALSE) +
    scale_fill_manual(values = info$colors, guide = "none") +
    scale_colour_manual(values = info$colors, guide = "none") +
    labs(title    = tp_d,
         subtitle = subtitle,
         x        = info$xlab,
         y        = paste(lipid_name, "intensity")) +
    theme_growell() +
    theme(plot.title    = element_text(size = FS_BASE, face = "bold"),
          plot.subtitle = element_text(size = FS_SMALL,
                                       colour = GROWELL_COLORS$grey_dark),
          axis.text.x   = element_text(size = FS_SMALL))
}

# ── 3-timepoint combined plot ─────────────────────────────
plot_lipid_3tp <- function(lipid_name, outcome, info) {
  plots <- lapply(seq_along(tp_names), function(i)
    plot_one_tp(lipid_name, i, outcome, info))

  wrap_plots(plots, ncol = ncols) +
    plot_annotation(
      title = paste0(lipid_name, " \u2014 ", info$xlab),
      theme = theme(plot.title = element_text(size = FS_TITLE, face = "bold",
                                               family = FONT_FAMILY))
    )
}

# ── Generate and save ─────────────────────────────────────
message("Generating lipid interest boxplots...")

for (outcome in names(lipids_of_interest)) {
  info   <- outcome_info[[outcome]]
  panels <- lipids_of_interest[[outcome]]

  for (panel_id in names(panels)) {
    lipid <- panels[[panel_id]]
    message(sprintf("  %s | %s | Panel %s: %s",
                    info$fig_id, outcome, panel_id, lipid))

    p <- plot_lipid_3tp(lipid, outcome, info)

    # Safe filename (replace spaces and colons)
    safe_name <- gsub("[ :]", "_", lipid)

    ggsave(file.path(out, paste0(info$fig_id, "_", panel_id, "_",
                                  safe_name, "_", plot_layout, ".pdf")),
           plot = p, width = panel_width, height = panel_height,
           units = "mm", device = cairo_pdf)

    save_panel(p, info$fig_id, panel_id,
               width_mm = panel_width, height_mm = panel_height)
  }
}

message("\nScript 11 complete → ", out)
message("Layout: ", plot_layout,
        " | To change: set plot_layout <- 'vertical' or 'horizontal'")