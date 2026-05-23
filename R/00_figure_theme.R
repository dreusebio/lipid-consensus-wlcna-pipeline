# R/00_figure_theme.R
# -------------------------------------------------------
# Global figure standards for npj Women's Health submission
#
# Requirements:
#   - Font: Arial (Helvetica fallback)
#   - DPI: 300 minimum (we use 300 for drafts, 600 for final)
#   - Font size: 8pt base (as recommended for print)
#   - Colorblind-safe palettes only
#   - No red/green combinations
#   - White background
#
# Usage: source("R/00_figure_theme.R") in any script
#        All plotting functions then use growell_* objects
# -------------------------------------------------------

# ── Font setup ────────────────────────────────────────────
# Register Arial/Helvetica for PDF and PNG output
library(ggplot2)
library(showtext)
library(sysfonts)

# Add Arial (falls back to Liberation Sans / Helvetica if not available)
tryCatch({
  font_add("Arial", regular = "Arial.ttf",
           bold = "Arial Bold.ttf",
           italic = "Arial Italic.ttf",
           bolditalic = "Arial Bold Italic.ttf")
  showtext_auto()
  FONT_FAMILY <- "Arial"
}, error = function(e) {
  tryCatch({
    FONT_FAMILY <<- "Helvetica"
  }, error = function(e2) {
    FONT_FAMILY <<- "Helvetica"
    message("Note: Arial not found, using Helvetica")
  })
})

# ── DPI settings ──────────────────────────────────────────
FIGURE_DPI        <- 300   # minimum per journal; use 600 for final submission
FIGURE_DPI_FINAL  <- 600

# ── Font sizes (8pt base as recommended for print figures) ──
FS_BASE   <- 8    # axis labels, legend text
FS_TITLE  <- 9    # panel titles
FS_SMALL  <- 6    # secondary labels, annotations
FS_LARGE  <- 10   # figure titles (avoid — use in legends only)

# ── Colorblind-safe palettes ──────────────────────────────
# Based on Wong (2011) Nature Methods colorblind-safe palette
# Avoids red/green combinations

GROWELL_COLORS <- list(

  # Primary two-group palette (APO, PPWR outcomes)
  # Royal blue = no adverse outcome, Scarlet = adverse outcome
  # Option 2: Scarlet & Royal Blue — colorblind safe, high contrast
  two_group = c(
    "No"  = "#1D3D8F",   # royal blue
    "Yes" = "#E63946"    # scarlet
  ),

  # Control vs Intervention
  treatment = c(
    "Control"      = "#1D3D8F",   # royal blue
    "Intervention" = "#E63946"    # scarlet
  ),

  # Three timepoints
  timepoints = c(
    "Baseline"   = "#009E73",   # teal/green
    "Wk36_38"    = "#CC79A7",   # mauve/pink
    "Postpartum" = "#F0E442"    # yellow (use with dark outline)
  ),
  timepoints_dark = c(
    "Baseline"   = "#005C42",
    "Wk36_38"    = "#8B4F70",
    "Postpartum" = "#A09810"
  ),

  # Lipid class palette (up to 10 classes)
  lipid_class = c(
    "Glycerophosphocholines" = "#0072B2",
    "Triacylglycerols"       = "#E69F00",
    "Sphingomyelins"         = "#009E73",
    "Ceramides"              = "#CC79A7",
    "Fatty acyls"            = "#56B4E9",
    "Lysophosphocholines"    = "#D55E00",
    "Phosphatidylethanolamines" = "#F0E442",
    "Phosphatidylinositols"  = "#000000",
    "Other"                  = "#999999",
    "Unknown"                = "#CCCCCC"
  ),
    # Diverging for module-trait correlation heatmaps
    # WGCNA-style blue-white-red palette
    # Blue  = negative correlation
    # White = zero/no correlation
    # Red   = positive correlation
    heatmap_low  = "blue",
    heatmap_mid  = "white",
    heatmap_high = "red",

  # Sequential (for enrichment, p-values)
  sequential = c("#EEF7FF","#BDD7EE","#7DB5D8","#3A82B8","#1C5B8A","#0C3357"),

  # Significance markers
  sig_color   = "#D55E00",   # vermillion for significant
  ns_color    = "#999999",   # grey for non-significant

  # Neutral greys
  grey_light  = "#F5F5F5",
  grey_mid    = "#CCCCCC",
  grey_dark   = "#666666"
)

# ── ggplot2 global theme ──────────────────────────────────
theme_growell <- function(base_size = FS_BASE,
                          base_family = FONT_FAMILY,
                          grid = "none") {

  t <- theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      # Panel
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.border      = element_rect(fill = NA, colour = "black",
                                       linewidth = 0.5),
      panel.grid.major  = if (grid == "both" || grid == "y")
                            element_line(colour = "#EEEEEE", linewidth = 0.3)
                          else element_blank(),
      panel.grid.minor  = element_blank(),

      # Axes
      axis.line         = element_blank(),
      axis.ticks        = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      axis.text         = element_text(size = FS_BASE, colour = "black",
                                       family = base_family),
      axis.title        = element_text(size = FS_BASE, colour = "black",
                                       family = base_family),

      # Legend
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key        = element_rect(fill = "white", colour = NA),
      legend.text       = element_text(size = FS_SMALL, family = base_family),
      legend.title      = element_text(size = FS_BASE, face = "bold",
                                       family = base_family),
      legend.key.size   = unit(8, "pt"),

      # Strip (facets)
      strip.background  = element_rect(fill = "#F0F0F0", colour = "black",
                                       linewidth = 0.5),
      strip.text        = element_text(size = FS_BASE, face = "bold",
                                       family = base_family),

      # Titles
      plot.title        = element_text(size = FS_TITLE, face = "bold",
                                       family = base_family, hjust = 0),
      plot.subtitle     = element_text(size = FS_BASE, family = base_family,
                                       hjust = 0),
      plot.background   = element_rect(fill = "white", colour = NA),
      plot.margin       = unit(c(3, 3, 3, 3), "pt")
    )
  t
}

# Set as ggplot2 default
theme_set(theme_growell())

# ── ColorRamp for heatmaps ────────────────────────────────
growell_heatmap_colors <- function(n = 100) {
  colorRampPalette(c(
    GROWELL_COLORS$heatmap_low,
    GROWELL_COLORS$heatmap_mid,
    GROWELL_COLORS$heatmap_high
  ))(n)
}

# ComplexHeatmap color function (colorblind-safe diverging)
growell_heatmap_col_fun <- function(range = c(-1, 0, 1)) {
  circlize::colorRamp2(
    range,
    c(GROWELL_COLORS$heatmap_low,
      GROWELL_COLORS$heatmap_mid,
      GROWELL_COLORS$heatmap_high)
  )
}

# ── Figure saving helper ──────────────────────────────────
# Saves individual panel files AND tracks them for later assembly
# Panel files are saved to results/figures/panels/<fig_id>/
# Final assembled figures go to results/figures/final/

FIGURE_REGISTRY <- list()  # tracks all panels in this session

save_panel <- function(plot_obj, fig_id, panel_id,
                       width_mm, height_mm,
                       dpi = FIGURE_DPI,
                       format = "pdf") {

  panel_dir <- file.path("results", "figures", "panels", fig_id)
  dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)

  filename <- paste0(fig_id, "_", panel_id, ".", format)
  filepath <- file.path(panel_dir, filename)

  if (format == "pdf") {
    # PDF for vector output (preferred by journal)
    ggsave(filepath, plot = plot_obj,
           width = width_mm, height = height_mm, units = "mm",
           device = cairo_pdf)
  } else {
    ggsave(filepath, plot = plot_obj,
           width = width_mm, height = height_mm, units = "mm",
           dpi = dpi)
  }

  # Also save PNG for easy viewing
  png_path <- file.path(panel_dir,
                        paste0(fig_id, "_", panel_id, ".png"))
  ggsave(png_path, plot = plot_obj,
         width = width_mm, height = height_mm, units = "mm",
         dpi = dpi)

  # Register panel for assembly
  FIGURE_REGISTRY[[paste0(fig_id, "_", panel_id)]] <<- list(
    fig_id   = fig_id,
    panel_id = panel_id,
    filepath = filepath,
    png_path = png_path,
    width_mm = width_mm,
    height_mm = height_mm
  )

  message(sprintf("  Saved panel %s_%s → %s", fig_id, panel_id, filepath))
  invisible(filepath)
}

# ── Figure dimensions (mm) for a two-column journal ───────
# npj Women's Health is online-only so size is flexible
# Standard Nature figure widths:
#   Single column: 89mm
#   1.5 column:    120mm
#   Full page:     183mm (max)
#   Max height:    247mm

FIG_WIDTH_SINGLE <- 89
FIG_WIDTH_HALF   <- 120
FIG_WIDTH_FULL   <- 183
FIG_HEIGHT_MAX   <- 247

message("Figure theme loaded: Helvetica, 300 DPI draft / 600 DPI final, WGCNA blue-white-red heatmap palette")
message("Use save_panel() to save individual panels for later assembly")

# # R/00_figure_theme.R
# # -------------------------------------------------------
# # Global figure standards for npj Women's Health submission
# #
# # Requirements:
# #   - Font: Arial (Helvetica fallback)
# #   - DPI: 300 minimum (we use 300 for drafts, 600 for final)
# #   - Font size: 8pt base (as recommended for print)
# #   - Colorblind-safe palettes only
# #   - No red/green combinations
# #   - White background
# #
# # Usage: source("R/00_figure_theme.R") in any script
# #        All plotting functions then use growell_* objects
# # -------------------------------------------------------

# # ── Font setup ────────────────────────────────────────────
# # Register Arial/Helvetica for PDF and PNG output
# library(ggplot2)
# library(showtext)
# library(sysfonts)

# # Add Arial (falls back to Liberation Sans / Helvetica if not available)
# tryCatch({
#   font_add("Arial", regular = "Arial.ttf",
#            bold = "Arial Bold.ttf",
#            italic = "Arial Italic.ttf",
#            bolditalic = "Arial Bold Italic.ttf")
#   showtext_auto()
#   FONT_FAMILY <- "Arial"
# }, error = function(e) {
#   tryCatch({
#     font_add_google("Arimo", "Arial")  # Google Fonts Arial-equivalent
#     showtext_auto()
#     FONT_FAMILY <- "Arial"
#   }, error = function(e2) {
#     FONT_FAMILY <<- "Helvetica"
#     message("Note: Arial not found, using Helvetica")
#   })
# })

# # ── DPI settings ──────────────────────────────────────────
# FIGURE_DPI        <- 300   # minimum per journal; use 600 for final submission
# FIGURE_DPI_FINAL  <- 600

# # ── Font sizes (8pt base as recommended for print figures) ──
# FS_BASE   <- 8    # axis labels, legend text
# FS_TITLE  <- 9    # panel titles
# FS_SMALL  <- 6    # secondary labels, annotations
# FS_LARGE  <- 10   # figure titles (avoid — use in legends only)

# # ── Colorblind-safe palettes ──────────────────────────────
# # Based on Wong (2011) Nature Methods colorblind-safe palette
# # Avoids red/green combinations

# GROWELL_COLORS <- list(

#   # Primary two-group palette (APO, PPWR outcomes)
#   # Royal blue = no adverse outcome, Scarlet = adverse outcome
#   # Option 2: Scarlet & Royal Blue — colorblind safe, high contrast
#   two_group = c(
#     "No"  = "#1D3D8F",   # royal blue
#     "Yes" = "#E63946"    # scarlet
#   ),

#   # Control vs Intervention
#   treatment = c(
#     "Control"      = "#1D3D8F",   # royal blue
#     "Intervention" = "#E63946"    # scarlet
#   ),

#   # Three timepoints
#   timepoints = c(
#     "Baseline"   = "#009E73",   # teal/green
#     "Wk36_38"    = "#CC79A7",   # mauve/pink
#     "Postpartum" = "#F0E442"    # yellow (use with dark outline)
#   ),
#   timepoints_dark = c(
#     "Baseline"   = "#005C42",
#     "Wk36_38"    = "#8B4F70",
#     "Postpartum" = "#A09810"
#   ),

#   # Lipid class palette (up to 10 classes)
#   lipid_class = c(
#     "Glycerophosphocholines" = "#0072B2",
#     "Triacylglycerols"       = "#E69F00",
#     "Sphingomyelins"         = "#009E73",
#     "Ceramides"              = "#CC79A7",
#     "Fatty acyls"            = "#56B4E9",
#     "Lysophosphocholines"    = "#D55E00",
#     "Phosphatidylethanolamines" = "#F0E442",
#     "Phosphatidylinositols"  = "#000000",
#     "Other"                  = "#999999",
#     "Unknown"                = "#CCCCCC"
#   ),
#     # Diverging for module-trait correlation heatmaps
#     # WGCNA-style blue-white-red palette
#     # Blue  = negative correlation
#     # White = zero/no correlation
#     # Red   = positive correlation
#     heatmap_low  = "blue",
#     heatmap_mid  = "white",
#     heatmap_high = "red",

#   # Sequential (for enrichment, p-values)
#   sequential = c("#EEF7FF","#BDD7EE","#7DB5D8","#3A82B8","#1C5B8A","#0C3357"),

#   # Significance markers
#   sig_color   = "#D55E00",   # vermillion for significant
#   ns_color    = "#999999",   # grey for non-significant

#   # Neutral greys
#   grey_light  = "#F5F5F5",
#   grey_mid    = "#CCCCCC",
#   grey_dark   = "#666666"
# )

# # ── ggplot2 global theme ──────────────────────────────────
# theme_growell <- function(base_size = FS_BASE,
#                           base_family = FONT_FAMILY,
#                           grid = "none") {

#   t <- theme_bw(base_size = base_size, base_family = base_family) +
#     theme(
#       # Panel
#       panel.background  = element_rect(fill = "white", colour = NA),
#       panel.border      = element_rect(fill = NA, colour = "black",
#                                        linewidth = 0.5),
#       panel.grid.major  = if (grid == "both" || grid == "y")
#                             element_line(colour = "#EEEEEE", linewidth = 0.3)
#                           else element_blank(),
#       panel.grid.minor  = element_blank(),

#       # Axes
#       axis.line         = element_blank(),
#       axis.ticks        = element_line(colour = "black", linewidth = 0.3),
#       axis.ticks.length = unit(2, "pt"),
#       axis.text         = element_text(size = FS_BASE, colour = "black",
#                                        family = base_family),
#       axis.title        = element_text(size = FS_BASE, colour = "black",
#                                        family = base_family),

#       # Legend
#       legend.background = element_rect(fill = "white", colour = NA),
#       legend.key        = element_rect(fill = "white", colour = NA),
#       legend.text       = element_text(size = FS_SMALL, family = base_family),
#       legend.title      = element_text(size = FS_BASE, face = "bold",
#                                        family = base_family),
#       legend.key.size   = unit(8, "pt"),

#       # Strip (facets)
#       strip.background  = element_rect(fill = "#F0F0F0", colour = "black",
#                                        linewidth = 0.5),
#       strip.text        = element_text(size = FS_BASE, face = "bold",
#                                        family = base_family),

#       # Titles
#       plot.title        = element_text(size = FS_TITLE, face = "bold",
#                                        family = base_family, hjust = 0),
#       plot.subtitle     = element_text(size = FS_BASE, family = base_family,
#                                        hjust = 0),
#       plot.background   = element_rect(fill = "white", colour = NA),
#       plot.margin       = unit(c(3, 3, 3, 3), "pt")
#     )
#   t
# }

# # Set as ggplot2 default
# theme_set(theme_growell())

# # ── ColorRamp for heatmaps ────────────────────────────────
# growell_heatmap_colors <- function(n = 100) {
#   colorRampPalette(c(
#     GROWELL_COLORS$heatmap_low,
#     GROWELL_COLORS$heatmap_mid,
#     GROWELL_COLORS$heatmap_high
#   ))(n)
# }

# # ComplexHeatmap color function (colorblind-safe diverging)
# growell_heatmap_col_fun <- function(range = c(-1, 0, 1)) {
#   circlize::colorRamp2(
#     range,
#     c(GROWELL_COLORS$heatmap_low,
#       GROWELL_COLORS$heatmap_mid,
#       GROWELL_COLORS$heatmap_high)
#   )
# }

# # ── Figure saving helper ──────────────────────────────────
# # Saves individual panel files AND tracks them for later assembly
# # Panel files are saved to results/figures/panels/<fig_id>/
# # Final assembled figures go to results/figures/final/

# FIGURE_REGISTRY <- list()  # tracks all panels in this session

# save_panel <- function(plot_obj, fig_id, panel_id,
#                        width_mm, height_mm,
#                        dpi = FIGURE_DPI,
#                        format = "pdf") {

#   panel_dir <- file.path("results", "figures", "panels", fig_id)
#   dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)

#   filename <- paste0(fig_id, "_", panel_id, ".", format)
#   filepath <- file.path(panel_dir, filename)

#   if (format == "pdf") {
#     # PDF for vector output (preferred by journal)
#     ggsave(filepath, plot = plot_obj,
#            width = width_mm, height = height_mm, units = "mm",
#            device = cairo_pdf)
#   } else {
#     ggsave(filepath, plot = plot_obj,
#            width = width_mm, height = height_mm, units = "mm",
#            dpi = dpi)
#   }

#   # Also save PNG for easy viewing
#   png_path <- file.path(panel_dir,
#                         paste0(fig_id, "_", panel_id, ".png"))
#   ggsave(png_path, plot = plot_obj,
#          width = width_mm, height = height_mm, units = "mm",
#          dpi = dpi)

#   # Register panel for assembly
#   FIGURE_REGISTRY[[paste0(fig_id, "_", panel_id)]] <<- list(
#     fig_id   = fig_id,
#     panel_id = panel_id,
#     filepath = filepath,
#     png_path = png_path,
#     width_mm = width_mm,
#     height_mm = height_mm
#   )

#   message(sprintf("  Saved panel %s_%s → %s", fig_id, panel_id, filepath))
#   invisible(filepath)
# }

# # ── Figure dimensions (mm) for a two-column journal ───────
# # npj Women's Health is online-only so size is flexible
# # Standard Nature figure widths:
# #   Single column: 89mm
# #   1.5 column:    120mm
# #   Full page:     183mm (max)
# #   Max height:    247mm

# FIG_WIDTH_SINGLE <- 89
# FIG_WIDTH_HALF   <- 120
# FIG_WIDTH_FULL   <- 183
# FIG_HEIGHT_MAX   <- 247

# message("Figure theme loaded: Arial/Helvetica, 300 DPI, colorblind-safe palettes")
# message("Use save_panel() to save individual panels for later assembly")