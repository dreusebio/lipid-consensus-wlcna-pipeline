# R/05_module_trait_correlations.R
# -------------------------------------------------------
# Module-trait correlations with configurable raw-p/FDR heatmaps
# Method: pearson/spearman from config.yml
#
# Outputs:
#   results/05_module_trait_cor/<run_tag>/<method>/
#   results/figures/panels/fig2/<heatmap_mode>/
#   results/figures/final/Figure2_direct_<trait_set>_<mode>.pdf/png/tiff
# -------------------------------------------------------

source("R/00_load_packages.R")
source("R/00_figure_theme.R")

# ── Helpers ──────────────────────────────────────────────
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

get_current_wgcna_run_tag <- function(cfg) {
  tag_file <- file.path(cfg$output$processed, "current_wgcna_run_tag.txt")
  
  if (length(tag_file) != 1 || is.na(tag_file) || tag_file == "") {
    stop("Invalid current_wgcna_run_tag.txt path. Check cfg$output$processed.")
  }
  
  if (!file.exists(tag_file)) {
    stop(
      "Missing current WGCNA run tag file: ", tag_file, "\n",
      "Run R/03_run_consensus_wgcna.R first."
    )
  }
  
  trimws(readLines(tag_file, warn = FALSE)[1])
}

cutoff_label <- function(x) {
  # Gives common folder names:
  # 0.05 -> 005
  # 0.10 -> 01
  if (abs(x - 0.05) < 1e-9) return("005")
  if (abs(x - 0.10) < 1e-9) return("01")
  out <- sprintf("%.3f", x)
  out <- sub("^0\\.", "", out)
  out <- sub("0+$", "", out)
  out
}

normalize_module_names <- function(x) {
  x <- as.character(x)
  x <- gsub("^ME", "", x)
  paste0("ME", x)
}

load_module_order <- function(cfg, wgcna_run_tag, current_modules) {
  
  module_order_file <- file.path(
    cfg$output$processed,
    paste0("module_order_by_hub_lipid_class_", wgcna_run_tag, ".rds")
  )
  
  if (!file.exists(module_order_file)) {
    message("No hub-lipid-class module order file found. Using current WGCNA module order.")
    return(current_modules)
  }
  
  message("Loading module order from: ", module_order_file)
  obj <- readRDS(module_order_file)
  
  if (is.data.frame(obj)) {
    if ("Module" %in% colnames(obj)) {
      module_order <- obj$Module
    } else if ("module" %in% colnames(obj)) {
      module_order <- obj$module
    } else {
      message("Module order file has no Module/module column. Using current order.")
      return(current_modules)
    }
  } else {
    module_order <- obj
  }
  
  module_order <- normalize_module_names(module_order)
  current_modules_norm <- normalize_module_names(current_modules)
  
  ordered_keep <- module_order[module_order %in% current_modules_norm]
  remaining <- current_modules_norm[!(current_modules_norm %in% ordered_keep)]
  
  final_order <- c(ordered_keep, remaining)
  final_order <- final_order[!duplicated(final_order)]
  
  message("  Ordered modules loaded: ", length(ordered_keep))
  message("  Remaining modules appended: ", length(remaining))
  
  final_order
}

# ── Load config and data ─────────────────────────────────
cfg <- yaml::read_yaml("config/config.yml")

method <- cfg$module_trait$method %||% "spearman"
method <- tolower(method)

wgcna_run_tag <- get_current_wgcna_run_tag(cfg)

message("WGCNA run: ", wgcna_run_tag, " | Method: ", toupper(method))

out <- file.path(cfg$output$s05, wgcna_run_tag, method)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

consensusMEs_ID <- readRDS(file.path(
  cfg$output$processed,
  paste0("consensusMEs_ID_", wgcna_run_tag, ".rds")
))

demographic_data <- readRDS(file.path(
  cfg$output$processed,
  "demographic_data_bmi.rds"
))

exprSize_ID <- readRDS(file.path(
  cfg$output$processed,
  "exprSize_ID.rds"
))

All_nSets_ID <- exprSize_ID$nSets

tp_labels  <- c("Baseline", "Wk36_38", "Postpartum")
tp_display <- c("10–16 weeks", "36–38 weeks", "Postpartum")
panel_ids  <- c("A", "B", "C")

# Keep only numeric phenotype traits for correlation
pheno_data <- demographic_data
pheno_data <- pheno_data[, sapply(pheno_data, is.numeric), drop = FALSE]

# Heatmap / FDR options
heatmap_modes <- unlist(cfg$module_trait$heatmap_significance %||% "raw")
heatmap_modes <- tolower(heatmap_modes)

raw_p_cutoff <- cfg$module_trait$raw_p_cutoff %||% 0.05
fdr_cutoff   <- cfg$module_trait$fdr_cutoff %||% 0.05
fdr_method   <- cfg$module_trait$fdr_method %||% "BH"
fdr_scope    <- cfg$module_trait$fdr_scope %||% "family"

message("Heatmap significance modes: ", paste(heatmap_modes, collapse = ", "))
message(
  "Raw p cutoff: ", raw_p_cutoff,
  " | FDR cutoff: ", fdr_cutoff,
  " | FDR method: ", fdr_method,
  " | FDR scope: ", fdr_scope
)

# ── Correlation function ─────────────────────────────────
compute_cor <- function(x, y, method) {
  
  if (method == "spearman") {
    
    cor_mat <- stats::cor(
      x,
      y,
      method = "spearman",
      use = "pairwise.complete.obs"
    )
    
    n <- nrow(x)
    
    # Approximate p-values using t approximation
    # Protect against r = +/- 1
    cor_for_test <- pmin(pmax(cor_mat, -0.999999), 0.999999)
    t_stat <- cor_for_test * sqrt((n - 2) / (1 - cor_for_test^2))
    pval_mat <- 2 * stats::pt(-abs(t_stat), df = n - 2)
    
  } else if (method == "pearson") {
    
    cor_mat <- WGCNA::cor(
      x,
      y,
      method = "pearson",
      use = "p"
    )
    
    pval_mat <- WGCNA::corPvalueFisher(
      cor_mat,
      nrow(x),
      twoSided = TRUE
    )
    
  } else {
    stop("Unsupported correlation method: ", method)
  }
  
  dimnames(pval_mat) <- dimnames(cor_mat)
  
  list(cor = cor_mat, pval = pval_mat)
}

# ── Trait families for FDR ───────────────────────────────
make_trait_family_vector <- function(traits, cfg) {
  
  families <- cfg$module_trait$fdr_families
  
  trait_family <- rep("Exploratory_other", length(traits))
  names(trait_family) <- traits
  
  if (!is.null(families)) {
    for (fam in names(families)) {
      fam_traits <- unlist(families[[fam]])
      fam_traits <- fam_traits[fam_traits %in% traits]
      trait_family[fam_traits] <- fam
    }
  }
  
  trait_family
}

apply_fdr <- function(p_mat, trait_family, method = "BH", scope = "family") {
  
  q_mat <- matrix(
    NA_real_,
    nrow = nrow(p_mat),
    ncol = ncol(p_mat),
    dimnames = dimnames(p_mat)
  )
  
  if (scope == "all_traits") {
    
    q_mat[] <- p.adjust(as.vector(p_mat), method = method)
    
  } else if (scope == "family") {
    
    for (fam in unique(trait_family)) {
      
      fam_traits <- names(trait_family)[trait_family == fam]
      fam_traits <- fam_traits[fam_traits %in% colnames(p_mat)]
      
      if (length(fam_traits) == 0) next
      
      p_sub <- p_mat[, fam_traits, drop = FALSE]
      
      q_sub <- matrix(
        p.adjust(as.vector(p_sub), method = method),
        nrow = nrow(p_sub),
        ncol = ncol(p_sub),
        dimnames = dimnames(p_sub)
      )
      
      q_mat[, fam_traits] <- q_sub
    }
    
  } else {
    stop("Unsupported fdr_scope: ", scope)
  }
  
  q_mat
}

# ── Heatmap function for all-trait panels ─────────────────
make_heatmap_panel <- function(cor_mat,
                               pval_mat,
                               qval_mat,
                               title_str,
                               pdf_path,
                               fig_panel_id,
                               heatmap_mode = c("raw", "fdr"),
                               module_order = NULL) {
  
  heatmap_mode <- match.arg(heatmap_mode)
  
  if (is.null(cor_mat) || nrow(cor_mat) == 0 || ncol(cor_mat) == 0) {
    return(invisible(NULL))
  }
  
  # Reorder modules, with TG-associated modules first if module_order exists
  if (!is.null(module_order)) {
    keep_order <- module_order[module_order %in% rownames(cor_mat)]
    remaining  <- rownames(cor_mat)[!(rownames(cor_mat) %in% keep_order)]
    final_order <- c(keep_order, remaining)
    
    cor_mat  <- cor_mat[final_order, , drop = FALSE]
    pval_mat <- pval_mat[final_order, , drop = FALSE]
    qval_mat <- qval_mat[final_order, , drop = FALSE]
  }
  
  # Significance matrix
  if (heatmap_mode == "fdr") {
    sig_mat <- qval_mat
    cutoff  <- fdr_cutoff
    sig_note <- paste0("BH q < ", cutoff)
  } else {
    sig_mat <- pval_mat
    cutoff  <- raw_p_cutoff
    sig_note <- paste0("p < ", cutoff)
  }
  
  textMatrix <- apply(sig_mat, 2, function(x) {
    sapply(x, function(v) {
      ifelse(!is.na(v) && v < cutoff, "*", "")
    })
  })
  
  dim(textMatrix) <- dim(cor_mat)
  
  row_labels <- gsub("^ME", "", rownames(cor_mat))
  
  # Dynamic panel size
  panel_width_mm  <- FIG_WIDTH_FULL
  panel_height_mm <- max(175, nrow(cor_mat) * 5.8)
  panel_height_mm <- min(panel_height_mm, 220)
  
  heatmap_cex_text <- 0.55
  heatmap_cex_y    <- 0.62
  heatmap_cex_x    <- 0.68
  
  # Save vector PDF
  grDevices::cairo_pdf(
    filename = pdf_path,
    width    = panel_width_mm / 25.4,
    height   = panel_height_mm / 25.4,
    family   = "Helvetica"
  )
  
  par(
    mar = c(8.5, 7.8, 1.8, 3.5),
    oma = c(0, 0, 0, 0),
    family = "Helvetica",
    xpd = NA
  )
  
  WGCNA::labeledHeatmap(
    Matrix        = cor_mat,
    xLabels       = colnames(cor_mat),
    yLabels       = rownames(cor_mat),
    ySymbols      = row_labels,
    colorLabels   = FALSE,
    colors        = WGCNA::blueWhiteRed(50),
    textMatrix    = textMatrix,
    setStdMargins = FALSE,
    cex.text      = heatmap_cex_text,
    textAdj       = c(0.5, 0.8),
    zlim          = c(-1, 1),
    main          = title_str,
    cex.lab.y     = heatmap_cex_y,
    cex.lab.x     = heatmap_cex_x,
    plotLegend    = TRUE,
    legendLabel   = ifelse(method == "pearson", "Pearson r", "Spearman \u03c1")
  )
  
  mtext(sig_note, side = 1, line = 7.2, adj = 1, cex = 0.65)
  
  dev.off()
  
  # Save 600 dpi PNG
  png_path <- gsub("\\.pdf$", ".png", pdf_path)
  
  png(
    filename = png_path,
    width    = panel_width_mm,
    height   = panel_height_mm,
    units    = "mm",
    res      = FIGURE_DPI_FINAL,
    type     = "cairo",
    family   = "Helvetica"
  )
  
  par(
    mar = c(8.5, 7.8, 1.8, 3.5),
    oma = c(0, 0, 0, 0),
    family = "Helvetica",
    xpd = NA
  )
  
  WGCNA::labeledHeatmap(
    Matrix        = cor_mat,
    xLabels       = colnames(cor_mat),
    yLabels       = rownames(cor_mat),
    ySymbols      = row_labels,
    colorLabels   = FALSE,
    colors        = WGCNA::blueWhiteRed(50),
    textMatrix    = textMatrix,
    setStdMargins = FALSE,
    cex.text      = heatmap_cex_text,
    textAdj       = c(0.5, 0.8),
    zlim          = c(-1, 1),
    main          = title_str,
    cex.lab.y     = heatmap_cex_y,
    cex.lab.x     = heatmap_cex_x,
    plotLegend    = TRUE,
    legendLabel   = ifelse(method == "pearson", "Pearson r", "Spearman \u03c1")
  )
  
  mtext(sig_note, side = 1, line = 7.2, adj = 1, cex = 0.65)
  
  dev.off()
  
  # Register panel in heatmap-specific panel folder
  heatmap_folder_name <- basename(dirname(pdf_path))
  panel_dir <- file.path("results", "figures", "panels", "fig2", heatmap_folder_name)
  dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)
  
  file.copy(
    png_path,
    file.path(panel_dir, paste0("fig2_", fig_panel_id, ".png")),
    overwrite = TRUE
  )
  
  file.copy(
    pdf_path,
    file.path(panel_dir, paste0("fig2_", fig_panel_id, ".pdf")),
    overwrite = TRUE
  )
  
  message(
    "  ", heatmap_mode, " panel fig2_", fig_panel_id,
    " saved | height: ", round(panel_height_mm, 1), " mm → ",
    dirname(pdf_path)
  )
  
  invisible(list(pdf = pdf_path, png = png_path))
}

# ── Direct journal-quality Figure 2: primary traits only ──
make_direct_figure2_heatmaps <- function(moduleTraitCor_ID,
                                         moduleTraitPvalue_ID,
                                         moduleTraitQvalue_ID,
                                         module_order = NULL,
                                         cfg,
                                         out_dir) {
  
  message("\nCreating direct journal-quality Figure 2 heatmaps...")
  
  trait_set    <- cfg$module_trait$figure2_trait_set %||% "primary"
  heatmap_mode <- cfg$module_trait$figure2_heatmap_mode %||% "fdr"
  
  if (trait_set == "primary") {
    traits_use <- unlist(cfg$module_trait$primary_figure2_traits)
  } else {
    traits_use <- colnames(moduleTraitCor_ID[[1]])
  }
  
  traits_use <- traits_use[traits_use %in% colnames(moduleTraitCor_ID[[1]])]
  
  if (length(traits_use) == 0) {
    warning("No direct Figure 2 traits found. Skipping direct Figure 2.")
    return(invisible(NULL))
  }
  
  message("  Figure 2 trait set: ", trait_set)
  message("  Traits used: ", paste(traits_use, collapse = ", "))
  message("  Significance mode: ", heatmap_mode)
  
  fig_dir <- file.path("results", "figures", "final")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  
  fig_pdf <- file.path(fig_dir, paste0("Figure2_direct_", trait_set, "_", heatmap_mode, ".pdf"))
  fig_png <- file.path(fig_dir, paste0("Figure2_direct_", trait_set, "_", heatmap_mode, ".png"))
  fig_tif <- file.path(fig_dir, paste0("Figure2_direct_", trait_set, "_", heatmap_mode, ".tiff"))
  
  panel_titles <- c("10–16 weeks", "36–38 weeks", "Postpartum")
  panel_letters <- c("A", "B", "C")
  
  n_modules <- nrow(moduleTraitCor_ID[[1]])
  
  fig_width_mm  <- FIG_WIDTH_FULL
  fig_height_mm <- FIG_HEIGHT_MAX
  
  draw_figure2 <- function() {
    
    layout(matrix(1:3, nrow = 3, ncol = 1))
    
    par(
      oma = c(0.4, 0.4, 0.4, 0.4),
      family = "Helvetica"
    )
    
    for (set in seq_along(moduleTraitCor_ID)) {
      
      cor_mat <- moduleTraitCor_ID[[set]]
      p_mat   <- moduleTraitPvalue_ID[[set]]
      q_mat   <- moduleTraitQvalue_ID[[set]]
      
      cor_mat <- cor_mat[, traits_use, drop = FALSE]
      p_mat   <- p_mat[, traits_use, drop = FALSE]
      q_mat   <- q_mat[, traits_use, drop = FALSE]
      
      # Module order
      if (!is.null(module_order)) {
        keep_order <- module_order[module_order %in% rownames(cor_mat)]
        remaining  <- rownames(cor_mat)[!(rownames(cor_mat) %in% keep_order)]
        final_order <- c(keep_order, remaining)
        
        cor_mat <- cor_mat[final_order, , drop = FALSE]
        p_mat   <- p_mat[final_order, , drop = FALSE]
        q_mat   <- q_mat[final_order, , drop = FALSE]
      }
      
      row_labels <- gsub("^ME", "", rownames(cor_mat))
      
      if (heatmap_mode == "fdr") {
        sig_mat <- q_mat
        cutoff <- fdr_cutoff
        sig_label <- paste0("* BH q < ", cutoff)
      } else {
        sig_mat <- p_mat
        cutoff <- raw_p_cutoff
        sig_label <- paste0("* p < ", cutoff)
      }
      
      textMatrix <- apply(sig_mat, 2, function(x) {
        sapply(x, function(v) {
          ifelse(!is.na(v) && v < cutoff, "*", "")
        })
      })
      
      dim(textMatrix) <- dim(cor_mat)
      
      par(
        mar = c(4.8, 7.4, 2.2, 4.2),
        family = "Helvetica"
      )
      
      WGCNA::labeledHeatmap(
        Matrix        = cor_mat,
        xLabels       = colnames(cor_mat),
        yLabels       = rownames(cor_mat),
        ySymbols      = row_labels,
        colorLabels   = FALSE,
        colors        = WGCNA::blueWhiteRed(50),
        textMatrix    = textMatrix,
        setStdMargins = FALSE,
        cex.text      = 0.78,
        textAdj       = c(0.5, 0.8),
        zlim          = c(-1, 1),
        main          = paste0(panel_titles[set], " — ", toupper(method)),
        cex.lab.y     = 0.62,
        cex.lab.x     = 0.82,
        plotLegend    = TRUE,
        legendLabel   = ifelse(method == "pearson", "Pearson r", "Spearman \u03c1")
      )
      
      mtext(
        panel_letters[set],
        side = 3,
        line = 0.45,
        adj = -0.08,
        cex = 1.2,
        font = 2
      )
      
      mtext(
        sig_label,
        side = 1,
        line = 4.0,
        adj = 1,
        cex = 0.65
      )
    }
  }
  

    convert_pdf_to_rasters <- function(pdf_file, png_file, tiff_file, dpi = 600) {
    
    message("Converting PDF to high-resolution PNG/TIFF...")
    
    # Check for pdftocairo
    pdftocairo_path <- Sys.which("pdftocairo")
    
    if (pdftocairo_path == "") {
        warning(
        "pdftocairo not found. Skipping PNG/TIFF conversion.\n",
        "Install poppler or use the PDF directly."
        )
        return(invisible(NULL))
    }
    
    # Temporary output prefix
    tmp_prefix <- tempfile("figure2_pdf_render_")
    
    # Render PDF to PNG at high DPI
    system2(
        pdftocairo_path,
        args = c(
        "-png",
        "-r", as.character(dpi),
        pdf_file,
        tmp_prefix
        )
    )
    
    rendered_png <- paste0(tmp_prefix, "-1.png")
    
    if (!file.exists(rendered_png)) {
        warning("PDF-to-PNG conversion failed. Expected file not found: ", rendered_png)
        return(invisible(NULL))
    }
    
    file.copy(rendered_png, png_file, overwrite = TRUE)
    
    # Convert PNG to TIFF with LZW compression using magick if available
    if (requireNamespace("magick", quietly = TRUE)) {
        img <- magick::image_read(rendered_png)
        magick::image_write(
        img,
        path = tiff_file,
        format = "tiff",
        compression = "lzw",
        density = paste0(dpi, "x", dpi)
        )
    } else {
        warning("magick not available. PNG created but TIFF skipped.")
    }
    
    message("  PNG:  ", png_file)
    message("  TIFF: ", tiff_file)
    }
  
  # Vector PDF
  grDevices::cairo_pdf(
    filename = fig_pdf,
    width = fig_width_mm / 25.4,
    height = fig_height_mm / 25.4,
    family = "Helvetica"
  )
  draw_figure2()
  dev.off()

    convert_pdf_to_rasters(
    pdf_file  = fig_pdf,
    png_file  = fig_png,
    tiff_file = fig_tif,
    dpi       = FIGURE_DPI_FINAL
    )
    # ── Convert vector PDF to high-resolution PNG/TIFF ───────
    # This preserves the good PDF layout and gives cleaner raster exports.

  message("Direct Figure 2 saved:")
  message("  PDF:  ", fig_pdf)
  message("  PNG:  ", fig_png)
  message("  TIFF: ", fig_tif)
  
  invisible(list(pdf = fig_pdf, png = fig_png, tiff = fig_tif))
}

# ── Compute correlations and FDR ──────────────────────────
message("\nComputing module-trait correlations (", method, ")...")

trait_family <- make_trait_family_vector(colnames(pheno_data), cfg)

message("Trait families:")
print(table(trait_family))

moduleTraitCor_ID    <- list()
moduleTraitPvalue_ID <- list()
moduleTraitQvalue_ID <- list()

# Load module order after seeing the module names
first_me <- consensusMEs_ID[[1]]$data
module_order <- load_module_order(
  cfg = cfg,
  wgcna_run_tag = wgcna_run_tag,
  current_modules = colnames(first_me)
)

for (set in 1:All_nSets_ID) {
  
  me <- consensusMEs_ID[[set]]$data
  tp <- tp_labels[set]
  
  ids <- intersect(rownames(me), rownames(pheno_data))
  message(tp, " | overlap: ", length(ids), " samples")
  
  if (length(ids) == 0) {
    stop("No overlapping sample IDs for ", tp)
  }
  
  me_sub <- me[ids, , drop = FALSE]
  pheno_sub <- pheno_data[ids, , drop = FALSE]
  
  res <- compute_cor(me_sub, pheno_sub, method)
  
  p_mat <- res$pval
  
  q_mat <- apply_fdr(
    p_mat = p_mat,
    trait_family = trait_family[colnames(p_mat)],
    method = fdr_method,
    scope = fdr_scope
  )
  
  moduleTraitCor_ID[[set]]    <- res$cor
  moduleTraitPvalue_ID[[set]] <- p_mat
  moduleTraitQvalue_ID[[set]] <- q_mat
  
  message("  FDR q<", fdr_cutoff, ": ", sum(q_mat < fdr_cutoff, na.rm = TRUE),
          " module-trait cells")
}

names(moduleTraitCor_ID)    <- tp_labels
names(moduleTraitPvalue_ID) <- tp_labels
names(moduleTraitQvalue_ID) <- tp_labels

# ── Save all-trait heatmaps and Excel files ───────────────
heatmap_folders <- character(0)

for (mode in heatmap_modes) {
  
  if (mode == "raw") {
    folder_name <- paste0("heatmaps_raw_p", cutoff_label(raw_p_cutoff))
  } else if (mode == "fdr") {
    folder_name <- paste0(
      "heatmaps_",
      fdr_scope,
      "_",
      fdr_method,
      "_q",
      cutoff_label(fdr_cutoff)
    )
  } else {
    stop("Unsupported heatmap_significance mode: ", mode)
  }
  
  heatmap_dir <- file.path(out, folder_name)
  dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)
  heatmap_folders <- c(heatmap_folders, heatmap_dir)
  
  for (set in 1:All_nSets_ID) {
    
    tp  <- tp_labels[set]
    tpd <- tp_display[set]
    
    make_heatmap_panel(
      cor_mat      = moduleTraitCor_ID[[set]],
      pval_mat     = moduleTraitPvalue_ID[[set]],
      qval_mat     = moduleTraitQvalue_ID[[set]],
      title_str    = paste0(tpd, " — ", toupper(method)),
      pdf_path     = file.path(heatmap_dir, paste0("heatmap_", tp, ".pdf")),
      fig_panel_id = panel_ids[set],
      heatmap_mode = mode,
      module_order = module_order
    )
  }
}

# Save Excel files with correlations, p-values, q-values
for (set in 1:All_nSets_ID) {
  
  tp <- tp_labels[set]
  
  wb <- openxlsx::createWorkbook()
  
  openxlsx::addWorksheet(wb, "Correlations")
  openxlsx::addWorksheet(wb, "P_values")
  openxlsx::addWorksheet(wb, "BH_Q_values")
  
  openxlsx::writeData(
    wb,
    "Correlations",
    round(moduleTraitCor_ID[[set]], 3),
    rowNames = TRUE
  )
  
  openxlsx::writeData(
    wb,
    "P_values",
    signif(moduleTraitPvalue_ID[[set]], 3),
    rowNames = TRUE
  )
  
  openxlsx::writeData(
    wb,
    "BH_Q_values",
    signif(moduleTraitQvalue_ID[[set]], 3),
    rowNames = TRUE
  )
  
  openxlsx::saveWorkbook(
    wb,
    file.path(out, paste0("ModuleTraitCor_", tp, ".xlsx")),
    overwrite = TRUE
  )
}

# ── Direct clean Figure 2 for main manuscript ─────────────
make_direct_figure2_heatmaps(
  moduleTraitCor_ID    = moduleTraitCor_ID,
  moduleTraitPvalue_ID = moduleTraitPvalue_ID,
  moduleTraitQvalue_ID = moduleTraitQvalue_ID,
  module_order         = module_order,
  cfg                  = cfg,
  out_dir              = out
)

# ── Save RDS ──────────────────────────────────────────────
rds_tag <- paste0(wgcna_run_tag, "_", method)

saveRDS(
  moduleTraitCor_ID,
  file.path(cfg$output$processed, paste0("moduleTraitCor_ID_", rds_tag, ".rds"))
)

saveRDS(
  moduleTraitPvalue_ID,
  file.path(cfg$output$processed, paste0("moduleTraitPvalue_ID_", rds_tag, ".rds"))
)

saveRDS(
  moduleTraitQvalue_ID,
  file.path(cfg$output$processed, paste0("moduleTraitQvalue_ID_", rds_tag, ".rds"))
)

writeLines(
  paste0(wgcna_run_tag, "|", method),
  file.path(cfg$output$processed, "current_analysis_tag.txt")
)

message("\nScript 05 complete → ", out)
message("Heatmap folders created:")
for (hf in unique(heatmap_folders)) {
  message("  - ", hf)
}
message("Figure 2 panels saved to results/figures/panels/fig2/<heatmap_mode>/")
message("Direct Figure 2 saved to results/figures/final/")
# # R/05_module_trait_correlations.R
# # -------------------------------------------------------
# # Unadjusted module-trait correlations (unadjusted only)
# # Method (pearson/spearman) from config.yml
# # Outputs → results/05_module_trait_cor/<run_tag>/<method>/
# # Figure 2 panels A, B, C → results/figures/panels/fig2/
# # -------------------------------------------------------
# source("R/00_load_packages.R")
# source("R/00_figure_theme.R")

# cfg    <- yaml::read_yaml("config/config.yml")
# method <- cfg$module_trait$method

# wgcna_run_tag <- trimws(readLines(
#   file.path(cfg$output$processed, "current_wgcna_run_tag.txt")))
# message("WGCNA run: ", wgcna_run_tag, " | Method: ", toupper(method))

# out <- file.path(cfg$output$s05, wgcna_run_tag, method)
# dir.create(out, showWarnings = FALSE, recursive = TRUE)

# consensusMEs_ID  <- readRDS(file.path(cfg$output$processed,
#   paste0("consensusMEs_ID_", wgcna_run_tag, ".rds")))
# demographic_data <- readRDS(file.path(cfg$output$processed, "demographic_data_bmi.rds"))
# exprSize_ID      <- readRDS(file.path(cfg$output$processed, "exprSize_ID.rds"))

# All_nSets_ID <- exprSize_ID$nSets
# tp_labels    <- c("Baseline", "Wk36_38", "Postpartum")
# tp_display   <- c("10–16 weeks", "36–38 weeks", "Postpartum")
# pheno_data   <- demographic_data

# # ── Correlation function ──────────────────────────────────
# compute_cor <- function(x, y, method) {
#   if (method == "spearman") {
#     cor_mat  <- cor(x, y, method = "spearman", use = "pairwise.complete.obs")
#     n        <- nrow(x)
#     t_stat   <- cor_mat * sqrt((n - 2) / (1 - cor_mat^2))
#     pval_mat <- 2 * pt(-abs(t_stat), df = n - 2)
#   } else {
#     cor_mat  <- WGCNA::cor(x, y, method = "pearson", use = "p")
#     pval_mat <- WGCNA::corPvalueFisher(cor_mat, nrow(x), twoSided = TRUE)
#   }
#   list(cor = cor_mat, pval = pval_mat)
# }

# # ── Heatmap function (publication-quality) ────────────────
# # Uses colorblind-safe blue-white-orange palette
# # Saves both to results/05 AND as figure panel
# make_heatmap_panel <- function(cor_mat, pval_mat, title_str,
#                                pdf_path, fig_panel_id) {

#   if (is.null(cor_mat) || nrow(cor_mat) == 0 || ncol(cor_mat) == 0)
#     return(invisible(NULL))

#   # Significance stars
#   textMatrix <- apply(pval_mat, 2, function(x)
#     sapply(x, function(p) ifelse(!is.na(p) && p < 0.05, "*", "")))
#   dim(textMatrix) <- dim(cor_mat)

#   # Clean module names (remove ME prefix)
#   row_labels <- gsub("^ME", "", rownames(cor_mat))

#   # ── Save to results/05 folder ──────────────────────────
#   pdf(pdf_path,
#       width  = FIG_WIDTH_FULL / 25.4,   # mm to inches
#       height = 120 / 25.4,
#       family = "Helvetica")
#   par(mar = c(12, 6, 2, 1), family = "Helvetica")
#   labeledHeatmap(
#     Matrix        = cor_mat,
#     xLabels       = colnames(cor_mat),
#     yLabels       = rownames(cor_mat),
#     ySymbols      = row_labels,
#     colorLabels   = FALSE,
#     colors        = blueWhiteRed(50),  # colorblind-safe
#     textMatrix    = textMatrix,
#     setStdMargins = FALSE,
#     cex.text      = FS_SMALL / 12,
#     textAdj       = c(0.5, 0.8),
#     zlim          = c(-1, 1),
#     main          = title_str,
#     cex.lab.y     = FS_BASE / 12,
#     cex.lab.x     = FS_BASE / 12,
#     plotLegend    = TRUE,
#     legendLabel   = ifelse(method == "pearson", "Pearson r", "Spearman \u03c1")
#   )
#   dev.off()

#   # ── Also save as figure panel PNG for assembly ─────────
#   png_path <- gsub("\\.pdf$", ".png", pdf_path)
#   png(png_path,
#       width  = FIG_WIDTH_FULL, height = 120, units = "mm",
#       res    = FIGURE_DPI, family = "Helvetica")
#   par(mar = c(12, 6, 2, 1), family = "Helvetica")
#   labeledHeatmap(
#     Matrix        = cor_mat,
#     xLabels       = colnames(cor_mat),
#     yLabels       = rownames(cor_mat),
#     ySymbols      = row_labels,
#     colorLabels   = FALSE,
#     colors        = growell_heatmap_colors(50),
#     textMatrix    = textMatrix,
#     setStdMargins = FALSE,
#     cex.text      = FS_SMALL / 12,
#     textAdj       = c(0.5, 0.8),
#     zlim          = c(-1, 1),
#     main          = title_str,
#     cex.lab.y     = FS_BASE / 12,
#     cex.lab.x     = FS_BASE / 12,
#     plotLegend    = TRUE,
#     legendLabel   = ifelse(method == "pearson", "Pearson r", "Spearman \u03c1")
#   )
#   dev.off()

#   # Register in figure panel directory
#   panel_dir <- file.path("results", "figures", "panels", "fig2")
#   dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)
#   file.copy(png_path, file.path(panel_dir, paste0("fig2_", fig_panel_id, ".png")),
#             overwrite = TRUE)
#   file.copy(pdf_path, file.path(panel_dir, paste0("fig2_", fig_panel_id, ".pdf")),
#             overwrite = TRUE)
#   message("  Panel fig2_", fig_panel_id, " saved")
# }

# # ── Compute correlations ──────────────────────────────────
# message("\nComputing module-trait correlations (", method, ")...")

# moduleTraitCor_ID    <- list()
# moduleTraitPvalue_ID <- list()

# for (set in 1:All_nSets_ID) {
#   me  <- consensusMEs_ID[[set]]$data
#   tp  <- tp_labels[set]
#   ids <- intersect(rownames(me), rownames(pheno_data))
#   message(tp, " | overlap: ", length(ids), " samples")
#   if (length(ids) == 0) stop("No overlapping sample IDs for ", tp)

#   res <- compute_cor(me[ids,, drop = FALSE], pheno_data[ids,, drop = FALSE], method)
#   moduleTraitCor_ID[[set]]    <- res$cor
#   moduleTraitPvalue_ID[[set]] <- res$pval
# }

# # ── Save heatmaps as Figure 2 panels A, B, C ─────────────
# panel_ids <- c("A", "B", "C")
# for (set in 1:All_nSets_ID) {
#   tp  <- tp_labels[set]
#   tpd <- tp_display[set]

#   make_heatmap_panel(
#     moduleTraitCor_ID[[set]],
#     moduleTraitPvalue_ID[[set]],
#     paste0(tpd, " — ", toupper(method)),
#     file.path(out, paste0("heatmap_", tp, ".pdf")),
#     panel_ids[set]
#   )

#   # Save Excel
#   wb <- createWorkbook()
#   addWorksheet(wb, "Correlations"); addWorksheet(wb, "P_values")
#   writeData(wb, "Correlations", round(moduleTraitCor_ID[[set]], 3), rowNames = TRUE)
#   writeData(wb, "P_values", signif(moduleTraitPvalue_ID[[set]], 3), rowNames = TRUE)
#   saveWorkbook(wb, file.path(out, paste0("ModuleTraitCor_", tp, ".xlsx")),
#                overwrite = TRUE)
# }

# # ── Save RDS ──────────────────────────────────────────────
# rds_tag <- paste0(wgcna_run_tag, "_", method)
# saveRDS(moduleTraitCor_ID,
#         file.path(cfg$output$processed,
#                   paste0("moduleTraitCor_ID_", rds_tag, ".rds")))
# saveRDS(moduleTraitPvalue_ID,
#         file.path(cfg$output$processed,
#                   paste0("moduleTraitPvalue_ID_", rds_tag, ".rds")))

# writeLines(paste0(wgcna_run_tag, "|", method),
#            file.path(cfg$output$processed, "current_analysis_tag.txt"))

# message("\nScript 05 complete → ", out)
# message("Figure 2 panels A/B/C saved to results/figures/panels/fig2/")