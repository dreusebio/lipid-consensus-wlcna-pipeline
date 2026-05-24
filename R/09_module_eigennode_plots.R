# R/09_module_eigennode_plots.R
# -------------------------------------------------------
# Module eigennode vs trait plots across 3 timepoints
# Figure 3 panels A-D → results/figures/panels/fig3/
#
# Layout controlled by: plot_layout
#   "vertical"   → 3 timepoints stacked top-to-bottom (ncol=1)
#                  89mm wide × 150mm tall per panel
#   "horizontal" → 3 timepoints left-to-right (nrow=1)
#                  183mm wide × 65mm tall per panel
# -------------------------------------------------------
source("R/00_load_packages.R")
source("R/00_figure_theme.R")

cfg <- yaml::read_yaml("config/config.yml")

# Read analysis tags
analysis_tag_raw <- trimws(readLines(
  file.path(cfg$output$processed, "current_analysis_tag.txt")))
analysis_parts <- strsplit(analysis_tag_raw, "\\|")[[1]]
wgcna_run_tag  <- analysis_parts[1]
method         <- analysis_parts[2]
rds_tag        <- paste0(wgcna_run_tag, "_", method)
message("WGCNA run: ", wgcna_run_tag, " | Method: ", method)

out <- file.path(cfg$output$s09, wgcna_run_tag, method)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ── Layout toggle ─────────────────────────────────────────
# "vertical"   → timepoints stacked top-to-bottom (89mm × 150mm)
# "horizontal" → timepoints left-to-right        (183mm × 65mm)
plot_layout  <- "horizontal"   # ← change to "vertical" if needed
panel_width  <- if (plot_layout == "vertical") FIG_WIDTH_SINGLE else FIG_WIDTH_FULL
panel_height <- if (plot_layout == "vertical") 150 else 65
panel_ncols  <- if (plot_layout == "vertical") 1 else 3
message("Layout: ", plot_layout)

consensusMEs_ID      <- readRDS(file.path(cfg$output$processed,
  paste0("consensusMEs_ID_", wgcna_run_tag, ".rds")))
moduleTraitCor_ID    <- readRDS(file.path(cfg$output$processed,
  paste0("moduleTraitCor_ID_", rds_tag, ".rds")))
moduleTraitPvalue_ID <- readRDS(file.path(cfg$output$processed,
  paste0("moduleTraitPvalue_ID_", rds_tag, ".rds")))
demographic_data     <- readRDS(file.path(cfg$output$processed, "demographic_data_bmi.rds"))

EigennodeDF <- as.data.frame(consensusMEs_ID)
traitData   <- list(
  Baseline   = demographic_data,
  TP36_38    = demographic_data,
  Postpartum = demographic_data
)

# ── Colorblind-safe two-group palette ─────────────────────
# Blue = no adverse outcome, Orange = adverse outcome (Wong 2011)
GROUP_COLORS_APO  <- c("No APO"  = "#1D3D8F",
                        "APO"     = "#E63946")
GROUP_COLORS_PPWR <- c("No PPWR" = "#1D3D8F",
                        "PPWR"    = "#E63946")
GROUP_COLORS_HDP  <- c("No HDP"  = "#1D3D8F",
                        "HDP"     = "#E63946")
GROUP_COLORS_GDM  <- c("No GDM"  = "#1D3D8F",
                        "GDM"     = "#E63946")
GROUP_COLORS_TRT  <- c("Control" = "#1D3D8F", "Intervention" = "#E63946")

# ── Single module-trait-timepoint plot ────────────────────
plot_me_trait_tp <- function(ME_vec, trait_vec, tp_label,
                              cor_val, pval,
                              trait_labels = NULL,
                              colors = NULL,
                              is_categorical = TRUE,
                              ylab = "Module eigennode",
                              xlab = "Group") {

  df <- na.omit(data.frame(Trait = trait_vec, ME = ME_vec))
  if (nrow(df) == 0) return(ggplot() + theme_void())

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
      levels(df$Trait) <- names(trait_labels)[match(levels(df$Trait),
                                                     as.character(trait_labels))]
    }
    p <- ggplot(df, aes(x = Trait, y = ME, fill = Trait)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7,
                   linewidth = 0.3, width = 0.5) +
      geom_jitter(aes(colour = Trait), width = 0.12,
                  size = 0.8, alpha = 0.6, show.legend = FALSE) +
      scale_fill_manual(values = colors, guide = "none") +
      scale_colour_manual(values = colors, guide = "none")
  } else {
    p <- ggplot(df, aes(x = Trait, y = ME)) +
      geom_point(colour = "#1D3D8F",
                 size = 1, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE,
                  colour = "#E63946",
                  linewidth = 0.5, fill = GROWELL_COLORS$grey_mid,
                  alpha = 0.2)
  }

  p + labs(title = tp_label, subtitle = subtitle,
           x = xlab, y = ylab) +
    theme_growell() +
    theme(plot.title    = element_text(size = FS_BASE, face = "bold"),
          plot.subtitle = element_text(size = FS_SMALL, colour = GROWELL_COLORS$grey_dark),
          axis.text.x   = element_text(size = FS_SMALL),
          legend.position = "none")
}

# ── Three-timepoint panel for one module × trait ──────────
plot_module_trait_3tp <- function(module, trait,
                                   trait_labels = NULL,
                                   colors = NULL,
                                   is_categorical = TRUE,
                                   ylab = NULL,
                                   xlab = trait) {

  tp_keys    <- c("Baseline", "TP36_38", "Postpartum")
  tp_display <- c("10\u201316 weeks", "36\u201338 weeks", "Postpartum")
  me_keys    <- paste0(c("Identified_Baseline",
                          "Identified_TP36_38weeks",
                          "Identified_Postpartum"),
                       ".data.ME", module)

  if (is.null(ylab)) ylab <- paste(module, "module eigennode")

  plots <- lapply(seq_along(tp_keys), function(i) {
    tp   <- tp_keys[i]
    me   <- EigennodeDF[[me_keys[i]]]
    tr   <- traitData[[tp]][[trait]]
    if (is.null(me) || is.null(tr)) return(ggplot() + theme_void())

    idx <- complete.cases(me, tr)
    me_sub <- me[idx]; tr_sub <- tr[idx]

    cor_val <- tryCatch(
      moduleTraitCor_ID[[i]][paste0("ME", module), trait],
      error = function(e) NA)
    pval <- tryCatch(
      moduleTraitPvalue_ID[[i]][paste0("ME", module), trait],
      error = function(e) NA)

    plot_me_trait_tp(me_sub, tr_sub, tp_display[i],
                     cor_val, pval, trait_labels, colors,
                     is_categorical, ylab, xlab)
  })

  wrap_plots(plots, ncol = 1) +
    plot_annotation(title = paste(module, "module \u2013", xlab),
                    theme = theme(
                      plot.title = element_text(size = FS_TITLE, face = "bold",
                                                family = FONT_FAMILY)))
}

# ── Trait info ────────────────────────────────────────────
traitInfo <- list(
  apo     = list(type="categorical", labels=c("No APO"=0,"APO"=1),
                 colors=GROUP_COLORS_APO,  xlab="APO status"),
  ppwr_e  = list(type="categorical", labels=c("No PPWR"=0,"PPWR"=1),
                 colors=GROUP_COLORS_PPWR, xlab="PPWR status"),
  apo_hdp = list(type="categorical", labels=c("No HDP"=0,"HDP"=1),
                 colors=GROUP_COLORS_HDP,  xlab="HDP status"),
  apo_gdm = list(type="categorical", labels=c("No GDM"=0,"GDM"=1),
                 colors=GROUP_COLORS_GDM,  xlab="GDM status"),
  group   = list(type="categorical", labels=c("Control"=0,"Intervention"=1),
                 colors=GROUP_COLORS_TRT,  xlab="Study arm"),
  bmi        = list(type="continuous", xlab="BMI"),
  bmi_base   = list(type="continuous", xlab="Baseline BMI"),
  gwg        = list(type="continuous", xlab="Gestational weight gain (kg)"),
  ppwr       = list(type="continuous", xlab="Postpartum weight retention (kg)"),
  weight_base= list(type="continuous", xlab="Baseline weight (kg)"),
  weight_wk36= list(type="continuous", xlab="Weight at 36\u201338 weeks (kg)"),
  weight_mth3= list(type="continuous", xlab="Weight at 3 months postpartum (kg)")
)

# ── Run all module × trait combinations ───────────────────
# Edit modules list based on your results from script 04
modules <- c("magenta","purple","brown","tan","grey60","red",
             "blue","lightcyan","lightyellow","cyan","black","midnightblue")
traits  <- names(traitInfo)

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
    fname <- file.path(out, paste0(module, "_vs_", trait, "_3TP.pdf"))
    ggsave(fname, plot = p,
           width = panel_width, height = panel_height, units = "mm",
           device = cairo_pdf)
    message("  ", module, " vs ", trait)
  }
}

# ── Save Figure 3 panels A-D ──────────────────────────────
# Edit these 4 combinations to be your key modules of interest
# after reviewing the heatmaps from script 05
# Currently set as placeholders — update once you know which modules matter

message("\nSaving Figure 3 panels (key modules of interest)...")
message("NOTE: Edit fig3_panels below after reviewing script 05 results")

fig3_panels <- list(
  A = list(module = modules[1], trait = "apo"),
  B = list(module = modules[1], trait = "ppwr_e"),
  C = list(module = modules[2], trait = "apo"),
  D = list(module = modules[2], trait = "ppwr_e")
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

  # Save panel
  save_panel(p, "fig3", panel_id,
             width_mm = panel_width, height_mm = panel_height)
  message("  Panel fig3_", panel_id, ": ",
          spec$module, " vs ", spec$trait)
}

message("\nScript 09 complete → ", out)
message("Layout: ", plot_layout, " | To change: set plot_layout <- 'vertical' or 'horizontal'")
message("Figure 3 panels saved to results/figures/panels/fig3/")
message("IMPORTANT: Update fig3_panels in this script with your key modules")
message("  after reviewing the heatmaps in results/05_module_trait_cor/")