# R/99_assemble_figures.R
# -------------------------------------------------------
# Assemble final publication figures from individual panels
#
# Figure layout:
#
# Figure 1 — Enrichment analysis (MetaboAnalyst)
#   A: Lipid class enrichment — all identified lipids
#   B: Pathway enrichment (if available)
#
# Figure 2 — WGCNA modules and module-trait correlations
#   A: Module-trait heatmap (Baseline)
#   B: Module-trait heatmap (36-38wk)
#   C: Module-trait heatmap (Postpartum)
#   D: Lipid class distribution by module
#
# Figure 3 — Module eigennode plots (key modules of interest)
#   A: Module vs outcome — timepoint 1
#   B: Module vs outcome — timepoint 2
#   C: Module vs outcome — timepoint 3
#   D: Module vs outcome — fourth panel (different module or outcome)
#
# Figure 4 — APO outcomes
#   A: Volcano/bar plot — differential lipids (APO)
#   B: PLS-DA biplot (APO)
#   C: Key lipid association 1 (APO)
#   D: Key lipid association 2 (APO)
#
# Figure 5 — PPWR outcomes
#   A: Volcano/bar plot — differential lipids (PPWR)
#   B: PLS-DA biplot (PPWR)
#   C: Key lipid association 1 (PPWR)
#   D: Key lipid association 2 (PPWR)
#
# Run after ALL analysis scripts have completed.
# -------------------------------------------------------

source("R/00_load_packages.R")
source("R/00_figure_theme.R")
library(patchwork)
library(png)
library(grid)
library(magick)

cfg           <- yaml::read_yaml("config/config.yml")
wgcna_run_tag <- trimws(readLines(
  file.path(cfg$output$processed, "current_wgcna_run_tag.txt")))
analysis_tag_raw <- trimws(readLines(
  file.path(cfg$output$processed, "current_analysis_tag.txt")))
analysis_parts <- strsplit(analysis_tag_raw, "\\|")[[1]]
wgcna_run_tag  <- analysis_parts[1]
method         <- analysis_parts[2]
rds_tag        <- paste0(wgcna_run_tag, "_", method)

final_dir <- file.path("results", "figures", "final")
panel_dir <- file.path("results", "figures", "panels")
dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)

message("Assembling figures for: ", rds_tag)

# ── Helper: load a saved panel PNG as a ggplot grob ──────
load_panel_png <- function(fig_id, panel_id) {
  path <- file.path(panel_dir, fig_id,
                    paste0(fig_id, "_", panel_id, ".png"))
  if (!file.exists(path)) {
    warning("Panel not found: ", path)
    return(ggplot() +
             annotate("text", x=0.5, y=0.5,
                      label=paste("Panel", panel_id, "not yet generated"),
                      size=3, colour="grey50") +
             theme_void())
  }
  img <- png::readPNG(path)
  ggplot() +
    annotation_raster(img, xmin=0, xmax=1, ymin=0, ymax=1) +
    coord_fixed(ratio = nrow(img)/ncol(img)) +
    theme_void()
}

# ── Helper: save assembled figure ────────────────────────
save_figure <- function(fig, fig_num, width_mm, height_mm,
                        dpi = FIGURE_DPI_FINAL) {
  tag  <- paste0("Figure", fig_num)
  pdf_path <- file.path(final_dir, paste0(tag, ".pdf"))
  png_path <- file.path(final_dir, paste0(tag, ".png"))

  ggsave(pdf_path, plot = fig,
         width = width_mm, height = height_mm, units = "mm",
         device = cairo_pdf)
  ggsave(png_path, plot = fig,
         width = width_mm, height = height_mm, units = "mm",
         dpi = dpi)

  message(sprintf("Saved %s: %s", tag, pdf_path))
  invisible(pdf_path)
}

# ── Figure 1: Enrichment analysis ────────────────────────
message("\nAssembling Figure 1 — Enrichment analysis...")

fig1_A <- load_panel_png("fig1", "A")  # lipid class enrichment
fig1_B <- load_panel_png("fig1", "B")  # pathway enrichment

fig1 <- (fig1_A | fig1_B) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(size = FS_TITLE, face = "bold",
                                           family = FONT_FAMILY))
  )

save_figure(fig1, 1, width_mm = FIG_WIDTH_FULL, height_mm = 100)

# ── Figure 2: WGCNA modules and module-trait heatmaps ────
# ── Figure 2: Module-trait heatmaps only: A, B, C ─────────
# ── Figure 2: Module-trait heatmaps A/B/C ────────────────
message("\nAssembling Figure 2 — full-width module-trait heatmaps A/B/C...")

fig2_heatmap_subdir <- if (!is.null(cfg$module_trait$figure2_heatmap_subdir)) {
  cfg$module_trait$figure2_heatmap_subdir
} else {
  "heatmaps_family_BH_q01"
}

fig2_panel_dir <- file.path(
  "results", "figures", "panels", "fig2", fig2_heatmap_subdir
)

fig2_paths <- file.path(
  fig2_panel_dir,
  paste0("fig2_", c("A", "B", "C"), ".png")
)

names(fig2_paths) <- c("A", "B", "C")

missing_panels <- fig2_paths[!file.exists(fig2_paths)]

if (length(missing_panels) > 0) {
  stop(
    "Missing Figure 2 panel PNG files:\n",
    paste(missing_panels, collapse = "\n"),
    "\nRun R/05_module_trait_correlations.R first."
  )
}

# Trim whitespace and standardize width
prepare_trimmed_panel <- function(path, panel_label, target_width_px = 5200) {
  
  message("  Preparing panel ", panel_label, ": ", path)
  
  img <- magick::image_read(path)
  
  # Trim white space around each heatmap
  img <- magick::image_trim(img, fuzz = 8)
  
  # Resize to common full width
  img <- magick::image_resize(img, paste0(target_width_px, "x"))
  
  # Save temporary trimmed PNG
  tmp_path <- tempfile(fileext = paste0("_fig2_", panel_label, ".png"))
  magick::image_write(img, tmp_path)
  
  tmp_path
}

trimmed_paths <- mapply(
  prepare_trimmed_panel,
  path = fig2_paths,
  panel_label = names(fig2_paths),
  SIMPLIFY = TRUE
)

# Read trimmed panels into R
panel_imgs <- lapply(trimmed_paths, png::readPNG)

# Get dimensions
panel_dims <- lapply(panel_imgs, function(x) {
  c(height = dim(x)[1], width = dim(x)[2])
})

target_width_px <- max(sapply(panel_dims, function(x) x["width"]))
panel_heights_px <- sapply(panel_dims, function(x) x["height"])

# Add a compact label strip before each heatmap
label_strip_px <- 170
gap_px <- 80

total_height_px <- sum(panel_heights_px) +
  length(panel_imgs) * label_strip_px +
  (length(panel_imgs) - 1) * gap_px

# Output paths
fig2_out_dir <- file.path("results", "figures", "final")
dir.create(fig2_out_dir, showWarnings = FALSE, recursive = TRUE)

fig2_png <- file.path(fig2_out_dir, paste0("Figure2_", fig2_heatmap_subdir, ".png"))
fig2_pdf <- file.path(fig2_out_dir, paste0("Figure2_", fig2_heatmap_subdir, ".pdf"))
fig2_tif <- file.path(fig2_out_dir, paste0("Figure2_", fig2_heatmap_subdir, ".tiff"))

# Helper to draw the final stacked figure using base R
draw_fig2 <- function() {
  
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(
    xlim = c(0, target_width_px),
    ylim = c(0, total_height_px),
    asp = 1
  )
  
  y_top <- total_height_px
  
  for (i in seq_along(panel_imgs)) {
    
    panel_label <- names(panel_imgs)[i]
    img <- panel_imgs[[i]]
    h <- dim(img)[1]
    w <- dim(img)[2]
    
    # Label strip
    y_top <- y_top - label_strip_px
    
    text(
      x = 40,
      y = y_top + label_strip_px * 0.55,
      labels = panel_label,
      adj = c(0, 0.5),
      cex = 3.8,
      font = 2,
      family = "sans"
    )
    
    # Center panel if needed
    x_left <- (target_width_px - w) / 2
    x_right <- x_left + w
    
    # Heatmap panel
    y_bottom <- y_top - h
    
    rasterImage(
      img,
      xleft = x_left,
      ybottom = y_bottom,
      xright = x_right,
      ytop = y_top,
      interpolate = TRUE
    )
    
    y_top <- y_bottom - gap_px
  }
}

# Save high-resolution PNG
png(
  filename = fig2_png,
  width = target_width_px,
  height = total_height_px,
  units = "px",
  res = 600,
  type = "cairo"
)
draw_fig2()
dev.off()

# Save PDF
grDevices::cairo_pdf(
  filename = fig2_pdf,
  width = target_width_px / 600,
  height = total_height_px / 600,
  family = "Helvetica"
)
draw_fig2()
dev.off()

# Save TIFF
tiff(
  filename = fig2_tif,
  width = target_width_px,
  height = total_height_px,
  units = "px",
  res = 600,
  compression = "lzw",
  type = "cairo"
)
draw_fig2()
dev.off()

message("Figure 2 saved:")
message("  PNG:  ", fig2_png)
message("  PDF:  ", fig2_pdf)
message("  TIFF: ", fig2_tif)

# ── Figure 3: Module eigennode plots ─────────────────────
message("\nAssembling Figure 3 — Module eigennode vs outcomes...")

fig3_A <- load_panel_png("fig3", "A")
fig3_B <- load_panel_png("fig3", "B")
fig3_C <- load_panel_png("fig3", "C")
fig3_D <- load_panel_png("fig3", "D")

fig3 <- (fig3_A | fig3_B) / (fig3_C | fig3_D) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(size = FS_TITLE, face = "bold",
                                           family = FONT_FAMILY))
  )

save_figure(fig3, 3, width_mm = FIG_WIDTH_FULL, height_mm = 160)

# ── Figure 4: APO outcomes ───────────────────────────────
message("\nAssembling Figure 4 — APO outcomes...")

fig4_A <- load_panel_png("fig4", "A")   # differential lipids
fig4_B <- load_panel_png("fig4", "B")   # PLS-DA biplot
fig4_C <- load_panel_png("fig4", "C")   # key lipid 1
fig4_D <- load_panel_png("fig4", "D")   # key lipid 2

fig4 <- (fig4_A | fig4_B) / (fig4_C | fig4_D) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(size = FS_TITLE, face = "bold",
                                           family = FONT_FAMILY))
  )

save_figure(fig4, 4, width_mm = FIG_WIDTH_FULL, height_mm = 160)

# ── Figure 5: PPWR outcomes ──────────────────────────────
message("\nAssembling Figure 5 — PPWR outcomes...")

fig5_A <- load_panel_png("fig5", "A")
fig5_B <- load_panel_png("fig5", "B")
fig5_C <- load_panel_png("fig5", "C")
fig5_D <- load_panel_png("fig5", "D")

fig5 <- (fig5_A | fig5_B) / (fig5_C | fig5_D) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(size = FS_TITLE, face = "bold",
                                           family = FONT_FAMILY))
  )

save_figure(fig5, 5, width_mm = FIG_WIDTH_FULL, height_mm = 160)

message("\n=== Figure assembly complete ===")
message("Final figures saved to: results/figures/final/")
message("To regenerate a single figure, rerun the relevant analysis script")
message("then rerun this assembly script.")