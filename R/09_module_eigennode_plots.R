# R/09_module_eigennode_plots.R
# -------------------------------------------------------
# Module eigennode vs trait plots across 3 timepoints
# Figure 3 panels A-D → results/figures/panels/fig3/
#
# Layout controlled by: plot_layout
#   "vertical"   → 3 timepoints stacked top-to-bottom
#                  85mm wide × 150mm tall per panel
#   "horizontal" → 3 timepoints left-to-right
#                  170mm wide × 55mm tall per panel
#
# To stack 4 panels in Illustrator at 230mm total height:
#   panel_height = 55mm (4 × 55mm + 3 × 3mm spacing = 229mm)
#
# To change which modules and traits are plotted, edit the "modules" and "traits" vectors in section 10 below.
# # ── 10. Generate all module × trait combinations ──────────
# # modules <- c("magenta", "purple", "brown", "tan", "grey60", "red",
# #              "blue", "lightcyan", "lightyellow", "cyan", "black",
# #              "midnightblue")
# # traits  <- names(traitInfo)

# # CHANGE TO ONLY THE MODULES YOU WANT:
# modules <- c("turquoise", "purple", "brown", "midnightblue")
# # add or remove any module names here

# # CHANGE TO ONLY THE TRAITS YOU WANT:
# traits <- c("apo", "ppwr_e")
# # add or remove any trait names here — must match keys in traitInfo

# -------------------------------------------------------

source("R/00_load_packages.R")
source("R/00_figure_theme.R")

library(ggplot2)
library(patchwork)

# ── 1. Helper: journal-compliant theme ───────────────────
# Defined at top so it is available to all plot functions.
# Helvetica 8pt throughout — meets npj Women's Health spec.
theme_eigennode_journal <- function() {
  theme_growell() +
  ggplot2::theme(
    # Panel border and grid
    panel.border      = ggplot2::element_rect(
                          colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.major  = ggplot2::element_line(
                          colour = "#EEEEEE", linewidth = 0.25),
    panel.grid.minor  = ggplot2::element_blank(),

    # Axes — 8pt bold
    axis.text.x       = ggplot2::element_text(
                          size = 8, face = "bold",
                          colour = "black", family = "Helvetica"),
    axis.text.y       = ggplot2::element_text(
                          size = 8, face = "bold",
                          colour = "black", family = "Helvetica"),
    axis.title.x      = ggplot2::element_text(
                          size = 8, face = "bold",
                          colour = "black", family = "Helvetica"),
    axis.title.y      = ggplot2::element_text(
                          size = 8, face = "bold",
                          colour = "black", family = "Helvetica"),
    axis.ticks        = ggplot2::element_line(
                          colour = "black", linewidth = 0.4),
    axis.ticks.length = ggplot2::unit(1.5, "mm"),

    # Title 8pt bold, subtitle 7pt plain
    plot.title        = ggplot2::element_text(
                          size = 8, face = "bold",
                          colour = "black", family = "Helvetica",
                          margin = ggplot2::margin(b = 2)),
    plot.subtitle     = ggplot2::element_text(
                          size = 7, face = "plain",
                          colour = "#444444", family = "Helvetica",
                          margin = ggplot2::margin(b = 3)),

    # Panel tag (a, b, c, d) — lower-case bold per journal spec
    plot.tag          = ggplot2::element_text(
                          size = 10, face = "bold",
                          colour = "black", family = "Helvetica"),
    plot.tag.position = c(0, 1),

    # Legend — 7pt, hidden by default for eigennode plots
    legend.text       = ggplot2::element_text(
                          size = 7, family = "Helvetica"),
    legend.title      = ggplot2::element_text(
                          size = 7, face = "bold", family = "Helvetica"),
    legend.key.size   = ggplot2::unit(3, "mm"),
    legend.position   = "none",

    # Tight margins to fill 170mm artboard
    plot.margin       = ggplot2::margin(t = 4, r = 4, b = 4, l = 4,
                                         unit = "mm")
  )
}

# ── 2. Helper: save vector PDF + 600 DPI TIFF ────────────
# cairo_pdf embeds Helvetica as a proper font object —
# text is fully selectable and editable in Adobe Illustrator.
# TIFF is provided as a backup raster for journal submission.
#feel free to adjust dimensions and resolution as needed for different figure sizes or journal requirements.
# save_eigennode_figure <- function(p, filepath_base,
#                                    width_mm  = 170,
#                                    height_mm = 55) {

save_eigennode_figure <- function(p, filepath_base,
                                   width_mm  = 170,
                                   height_mm = 55) {

  # PDF — Illustrator-editable vector
  grDevices::cairo_pdf(
    filename = paste0(filepath_base, ".pdf"),
    width    = width_mm  / 25.4,
    height   = height_mm / 25.4,
    family   = "Helvetica",
    onefile  = FALSE,
    bg       = "white"
  )
  print(p)
  grDevices::dev.off()

  # TIFF — 600 DPI LZW raster backup
  grDevices::tiff(
    filename    = paste0(filepath_base, ".tiff"),
    width       = width_mm,
    height      = height_mm,
    units       = "mm",
    res         = 600,
    compression = "lzw",
    type        = "cairo",
    family      = "Helvetica",
    bg          = "white"
  )
  print(p)
  grDevices::dev.off()

  message("    PDF + TIFF saved: ", basename(filepath_base))
}

# ── 3. Config and tags ────────────────────────────────────
cfg <- yaml::read_yaml("config/config.yml")

analysis_tag_raw <- trimws(readLines(
  file.path(cfg$output$processed, "current_analysis_tag.txt")))
analysis_parts <- strsplit(analysis_tag_raw, "\\|")[[1]]
wgcna_run_tag  <- analysis_parts[1]
method         <- tolower(trimws(analysis_parts[2]))
rds_tag        <- paste0(wgcna_run_tag, "_", method)
message("WGCNA run: ", wgcna_run_tag, " | Method: ", method)

out <- file.path(cfg$output$s09, wgcna_run_tag, method)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ── 4. Layout toggle ─────────────────────────────────────
# "horizontal" → 3 timepoints left-to-right, 170mm × 55mm
# "vertical"   → 3 timepoints stacked,       85mm  × 150mm
#
# For 4 panels stacked in Illustrator at 230mm total height:
#   panel_height = 55mm  (4 × 55 + 3 × 3mm gap = 229mm ≈ 230mm)
#
plot_layout  <- "horizontal"   # ← "horizontal" or "vertical"

panel_width  <- if (plot_layout == "horizontal") FIG_WIDTH_FULL   else FIG_WIDTH_SINGLE
panel_height <- if (plot_layout == "horizontal") 55               else 150
panel_ncols  <- if (plot_layout == "horizontal") 3                else 1

message("Layout: ", plot_layout,
        " | Size: ", panel_width, "mm × ", panel_height, "mm")

# ── 5. Load data ──────────────────────────────────────────
consensusMEs_ID      <- readRDS(file.path(cfg$output$processed,
  paste0("consensusMEs_ID_", wgcna_run_tag, ".rds")))
moduleTraitCor_ID    <- readRDS(file.path(cfg$output$processed,
  paste0("moduleTraitCor_ID_", rds_tag, ".rds")))
moduleTraitPvalue_ID <- readRDS(file.path(cfg$output$processed,
  paste0("moduleTraitPvalue_ID_", rds_tag, ".rds")))
demographic_data     <- readRDS(file.path(cfg$output$processed,
  "demographic_data_bmi.rds"))

EigennodeDF <- as.data.frame(consensusMEs_ID)
traitData   <- list(
  Baseline   = demographic_data,
  TP36_38    = demographic_data,
  Postpartum = demographic_data
)

# ── 6. Colour palettes ────────────────────────────────────
# Colorblind-safe: blue = no adverse outcome, red = adverse outcome
GROUP_COLORS_APO  <- c("No APO"       = "#1D3D8F", "APO"          = "#E63946")
GROUP_COLORS_PPWR <- c("No PPWR"      = "#1D3D8F", "PPWR"         = "#E63946")
GROUP_COLORS_HDP  <- c("No HDP"       = "#1D3D8F", "HDP"          = "#E63946")
GROUP_COLORS_GDM  <- c("No GDM"       = "#1D3D8F", "GDM"          = "#E63946")
GROUP_COLORS_TRT  <- c("Control"      = "#1D3D8F", "Intervention"  = "#E63946")

# ── 7. Single module-trait-timepoint plot ─────────────────
plot_me_trait_tp <- function(ME_vec, trait_vec, tp_label,
                              cor_val, pval,
                              trait_labels   = NULL,
                              colors         = NULL,
                              is_categorical = TRUE,
                              ylab           = "Module eigennode",
                              xlab           = "Group") {

  df <- na.omit(data.frame(Trait = trait_vec, ME = ME_vec))
  if (nrow(df) == 0) return(ggplot2::ggplot() + ggplot2::theme_void())

  stars <- ifelse(pval < 0.001, "***",
           ifelse(pval < 0.01,  "**",
           ifelse(pval < 0.05,  "*", "ns")))

  subtitle <- if (is_categorical)
    paste0("p = ", signif(pval, 2), " ", stars)
  else
    paste0(toupper(method), " \u03c1 = ", signif(cor_val, 2),
           ", p = ", signif(pval, 2), " ", stars)

  if (is_categorical) {
    df$Trait <- factor(df$Trait)
    if (!is.null(trait_labels)) {
      levels(df$Trait) <- names(trait_labels)[
        match(levels(df$Trait), as.character(trait_labels))]
    }
    p <- ggplot2::ggplot(df,
           ggplot2::aes(x = Trait, y = ME, fill = Trait)) +
      ggplot2::geom_boxplot(
        outlier.shape = NA, alpha = 0.7,
        linewidth = 0.3, width = 0.5) +
      ggplot2::geom_jitter(
        ggplot2::aes(colour = Trait),
        width = 0.12, size = 0.8, alpha = 0.6,
        show.legend = FALSE) +
      ggplot2::scale_fill_manual(  values = colors, guide = "none") +
      ggplot2::scale_colour_manual(values = colors, guide = "none")

  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = Trait, y = ME)) +
      ggplot2::geom_point(colour = "#1D3D8F", size = 1, alpha = 0.7) +
      ggplot2::geom_smooth(
        method = "lm", se = TRUE,
        colour = "#E63946", linewidth = 0.5,
        fill   = GROWELL_COLORS$grey_mid, alpha = 0.2)
  }

  # Apply theme — theme_eigennode_journal() called here,
  # NOT defined here
  p +
    ggplot2::labs(
      title    = tp_label,
      subtitle = subtitle,
      x        = xlab,
      y        = ylab) +
    theme_eigennode_journal()
}

# ── 8. Three-timepoint panel for one module × trait ───────
plot_module_trait_3tp <- function(module, trait,
                                   trait_labels   = NULL,
                                   colors         = NULL,
                                   is_categorical = TRUE,
                                   ylab           = NULL,
                                   xlab           = trait) {

  tp_keys    <- c("Baseline", "TP36_38", "Postpartum")
  tp_display <- c("10-16 weeks", "36-38 weeks", "Postpartum")
  me_keys    <- paste0(
    c("Identified_Baseline",
      "Identified_TP36_38weeks",
      "Identified_Postpartum"),
    ".data.ME", module)

  if (is.null(ylab)) ylab <- paste(module, "module eigennode")

  plots <- lapply(seq_along(tp_keys), function(i) {
    tp  <- tp_keys[i]
    me  <- EigennodeDF[[me_keys[i]]]
    tr  <- traitData[[tp]][[trait]]
    if (is.null(me) || is.null(tr))
      return(ggplot2::ggplot() + ggplot2::theme_void())

    idx    <- complete.cases(me, tr)
    me_sub <- me[idx]
    tr_sub <- tr[idx]

    cor_val <- tryCatch(
      moduleTraitCor_ID[[i]][paste0("ME", module), trait],
      error = function(e) NA_real_)
    pval <- tryCatch(
      moduleTraitPvalue_ID[[i]][paste0("ME", module), trait],
      error = function(e) NA_real_)

    plot_me_trait_tp(
      ME_vec         = me_sub,
      trait_vec      = tr_sub,
      tp_label       = tp_display[i],
      cor_val        = cor_val,
      pval           = pval,
      trait_labels   = trait_labels,
      colors         = colors,
      is_categorical = is_categorical,
      ylab           = ylab,
      xlab           = xlab)
  })

  patchwork::wrap_plots(plots, ncol = panel_ncols) +
    patchwork::plot_annotation(
      title = paste(module, "module -", xlab),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size   = 9, face   = "bold",
          family = "Helvetica",
          margin = ggplot2::margin(b = 3))
      )
    )
}

# ── 9. Trait definitions ──────────────────────────────────
traitInfo <- list(
  apo      = list(type = "categorical",
                  labels = c("No APO" = 0, "APO" = 1),
                  colors = GROUP_COLORS_APO,
                  xlab   = "APO status"),
  ppwr_e   = list(type = "categorical",
                  labels = c("No PPWR" = 0, "PPWR" = 1),
                  colors = GROUP_COLORS_PPWR,
                  xlab   = "PPWR status"),
  apo_hdp  = list(type = "categorical",
                  labels = c("No HDP" = 0, "HDP" = 1),
                  colors = GROUP_COLORS_HDP,
                  xlab   = "HDP status"),
  apo_gdm  = list(type = "categorical",
                  labels = c("No GDM" = 0, "GDM" = 1),
                  colors = GROUP_COLORS_GDM,
                  xlab   = "GDM status"),
  group    = list(type = "categorical",
                  labels = c("Control" = 0, "Intervention" = 1),
                  colors = GROUP_COLORS_TRT,
                  xlab   = "Study arm"),
  bmi         = list(type = "continuous", xlab = "BMI"),
  bmi_base    = list(type = "continuous", xlab = "Baseline BMI"),
  gwg         = list(type = "continuous",
                     xlab = "Gestational weight gain (kg)"),
  ppwr        = list(type = "continuous",
                     xlab = "Postpartum weight retention (kg)"),
  weight_base = list(type = "continuous",
                     xlab = "Baseline weight (kg)"),
  weight_wk36 = list(type = "continuous",
                     xlab = "Weight at 36\u201338 weeks (kg)"),
  weight_mth3 = list(type = "continuous",
                     xlab = "Weight at 3 months postpartum (kg)")
)

# ── 10. Generate all module × trait combinations ──────────
# modules <- c("magenta", "purple", "brown", "tan", "grey60", "red",
#              "blue", "lightcyan", "lightyellow", "cyan", "black",
#              "midnightblue")
# traits  <- names(traitInfo)

# CHANGE TO ONLY THE MODULES YOU WANT:
modules <- c("turquoise", "purple","grey60", "darkorange", "lightcyan")
# add or remove any module names here

# CHANGE TO ONLY THE TRAITS YOU WANT:
traits <- c("apo", "ppwr_e","ppwr")
# add or remove any trait names here — must match keys in traitInfo

message("Generating module-trait plots...")

for (module in modules) {
  for (trait in traits) {
    info   <- traitInfo[[trait]]
    is_cat <- info$type == "categorical"

    p <- plot_module_trait_3tp(
      module         = module,
      trait          = trait,
      trait_labels   = if (is_cat) info$labels else NULL,
      colors         = if (is_cat) info$colors else NULL,
      is_categorical = is_cat,
      xlab           = info$xlab
    )

    # save_eigennode_figure() called here — defined at top of script
    fname_base <- file.path(out,
                             paste0(module, "_vs_", trait, "_3TP"))
    save_eigennode_figure(p, fname_base,
                           width_mm  = panel_width,
                           height_mm = panel_height)
  }
}

# ── 11. Save Figure 3 key panels ─────────────────────────
# Edit fig3_panels to match your key modules after reviewing
# the heatmaps from script 08. Currently set as placeholders.
# Panel labels are lower-case bold (a, b, c, d) per journal spec.

message("\nSaving Figure 3 panels (key modules of interest)...")
message("NOTE: Update fig3_panels below after reviewing script 08 results")

fig3_panels <- list(
  a = list(module = modules[1], trait = "apo"),
  b = list(module = modules[1], trait = "ppwr_e"),
  c = list(module = modules[2], trait = "apo"),
  d = list(module = modules[2], trait = "ppwr_e")
)

panel_dir <- file.path("results", "figures", "panels", "fig3")
dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)

for (panel_id in names(fig3_panels)) {
  spec   <- fig3_panels[[panel_id]]
  info   <- traitInfo[[spec$trait]]
  is_cat <- info$type == "categorical"

  p <- plot_module_trait_3tp(
    module         = spec$module,
    trait          = spec$trait,
    trait_labels   = if (is_cat) info$labels else NULL,
    colors         = if (is_cat) info$colors else NULL,
    is_categorical = is_cat,
    xlab           = info$xlab
  )

  # Add lower-case bold panel label per journal spec
  p_labeled <- p +
    patchwork::plot_annotation(
      tag_levels = list(c(panel_id)),   # already lower-case (a, b, c, d)
      theme = ggplot2::theme(
        plot.tag = ggplot2::element_text(
          size   = 10,
          face   = "bold",
          family = "Helvetica"
        )
      )
    )

  panel_base <- file.path(panel_dir,
                           paste0("fig3_", panel_id))
  save_eigennode_figure(p_labeled, panel_base,
                         width_mm  = panel_width,
                         height_mm = panel_height)

  message("  Panel fig3_", panel_id, ": ",
          spec$module, " vs ", spec$trait)
}

message("\nScript 09 complete → ", out)
message("Layout: ", plot_layout,
        " | Panel size: ", panel_width, "mm × ", panel_height, "mm")
message("Figures saved as PDF (Illustrator-editable) + TIFF (600 DPI)")
message("Figure 3 panels: results/figures/panels/fig3/")
message("IMPORTANT: Update fig3_panels above with your key modules")
message("  after reviewing heatmaps in results/08_module_trait_cor/")