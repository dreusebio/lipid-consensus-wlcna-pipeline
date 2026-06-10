# R/08_module_trait_correlations.R
# -------------------------------------------------------
# Module-trait correlations with configurable raw-p/FDR heatmaps
# Adds lipid-class block strip beside module rows for direct Figure 2
# Orders TG-associated modules first
# Allows lower-case or upper-case panel labels from config.yml
# Method: pearson/spearman from config.yml
#
# Outputs:
#   results/05_module_trait_cor/<run_tag>/<method>/
#   results/figures/panels/fig2/<heatmap_mode>/
#   results/figures/final/Figure2_direct_<trait_set>_<mode>.pdf/png/tiff
# -------------------------------------------------------

source("R/00_load_packages.R")
source("R/00_figure_theme.R")

library(ggplot2)
library(grid)

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

# ── Lipid-class colors ───────────────────────────────────
# These are lipid-class colors, not WGCNA module colors.
# Grey is reserved only for the WGCNA grey module.

LIPID_CLASS_COLORS <- c(
  "TG"          = "#E63946",  # triglycerides / triacylglycerols
  "DG"          = "#F4A261",  # diacylglycerols
  "CE"          = "#8E44AD",  # cholesteryl esters
  "PC/LPC"      = "#1D3D8F",  # phosphatidylcholines / lysophosphatidylcholines
  "PE/LPE"      = "#F4A261",  # phosphatidylethanolamines / lysophosphatidylethanolamines
  "SM"          = "#2A9D8F",  # sphingomyelins
  "Cer"         = "#E9C46A",  # ceramides
  "FA"          = "#795548",  # fatty acyls
  "Oth"       = "#FFFFFF",  # other lipid class; not grey
  "Unknown"     = "#2a756cfb",  # unclassified
  "Grey module" = "#874788fb"   # WGCNA grey module only
)

LIPID_CLASS_PRIORITY <- c(
  "TG",
  "DG",
  "CE",
  "PC/LPC",
  "PE/LPE",
  "SM",
  "Cer",
  "FA",
  "Oth",
  "Unknown",
  "Grey module"
)

simplify_lipid_class <- function(x) {
  x <- as.character(x)
  
  # Handles labels such as:
  # "01_TG-associated", "05_PC/LPC", "06_Sphingomyelins", etc.
  x_clean <- x
  x_clean <- gsub("^[0-9]+_", "", x_clean)
  x_clean <- gsub("-associated", "", x_clean, ignore.case = TRUE)
  x_clean <- trimws(x_clean)
  
  out <- rep("Unknown", length(x_clean))
  
  out[grepl("^TG$|triacyl|triglycer", x_clean, ignore.case = TRUE)] <- "TG"
  out[grepl("^DG$|diacyl", x_clean, ignore.case = TRUE)] <- "DG"
  out[grepl("^CE$|cholesteryl|cholesterol ester", x_clean, ignore.case = TRUE)] <- "CE"
  out[grepl("PC/LPC|phosphatidylcholines|phosphocholine|\\bPC\\b|\\bLPC\\b",
            x_clean, ignore.case = TRUE)] <- "PC/LPC"
  out[grepl("PE/LPE|phosphatidylethanolamines|phosphoethanolamine|\\bPE\\b|\\bLPE\\b",
            x_clean, ignore.case = TRUE)] <- "PE/LPE"
  out[grepl("^SM$|sphingomyelin", x_clean, ignore.case = TRUE)] <- "SM"
  out[grepl("^Cer$|ceramide", x_clean, ignore.case = TRUE)] <- "Cer"
  out[grepl("fatty-acyl|fatty acyl|fatty", x_clean, ignore.case = TRUE)] <- "FA"
  out[grepl("mixed", x_clean, ignore.case = TRUE)] <- "Oth"
  out[grepl("other", x_clean, ignore.case = TRUE)] <- "Oth"
  
  out
}

# ── Module order and lipid-class annotation ───────────────

load_module_order_file <- function(cfg, wgcna_run_tag, current_modules) {
  
  module_order_file <- file.path(
    cfg$output$processed,
    paste0("module_order_by_hub_lipid_class_", wgcna_run_tag, ".rds")
  )
  
  current_modules_norm <- normalize_module_names(current_modules)
  
  if (!file.exists(module_order_file)) {
    message("No hub-lipid-class module order file found. Using current WGCNA module order.")
    return(current_modules_norm)
  }
  
  message("Loading module order from: ", module_order_file)
  obj <- readRDS(module_order_file)
  
  if (is.data.frame(obj)) {
    
    if (!("Module" %in% colnames(obj))) {
      stop("Module order file does not contain a 'Module' column.")
    }
    
    if ("PlotOrder" %in% colnames(obj)) {
      obj <- obj[order(obj$PlotOrder), , drop = FALSE]
    }
    
    module_order <- obj$Module
    
  } else {
    module_order <- obj
  }
  
  module_order <- normalize_module_names(module_order)
  
  ordered_keep <- module_order[module_order %in% current_modules_norm]
  remaining <- current_modules_norm[!(current_modules_norm %in% ordered_keep)]
  
  final_order <- c(ordered_keep, remaining)
  final_order <- final_order[!duplicated(final_order)]
  
  message("  Ordered modules loaded: ", length(ordered_keep))
  message("  Remaining modules appended: ", length(remaining))
  
  final_order
}

load_module_class_annotation <- function(cfg, wgcna_run_tag, current_modules) {
  
  module_class_file <- file.path(
    cfg$output$processed,
    paste0("module_order_by_hub_lipid_class_", wgcna_run_tag, ".rds")
  )
  
  current_modules_norm <- normalize_module_names(current_modules)
  
  default_df <- data.frame(
    Module     = current_modules_norm,
    LipidClass = "Unknown",
    ClassGroup = "Unknown",
    ClassColor = unname(LIPID_CLASS_COLORS["Unknown"]),
    stringsAsFactors = FALSE
  )
  
  if (!file.exists(module_class_file)) {
    message("No module lipid-class annotation file found. Using white unclassified strip.")
    return(default_df)
  }
  
  message("Loading module lipid-class annotation from: ", module_class_file)
  obj <- readRDS(module_class_file)
  
  if (!is.data.frame(obj)) {
    message("Module lipid-class annotation is not a data frame. Using white unclassified strip.")
    return(default_df)
  }
  
  if (!("Module" %in% colnames(obj))) {
    message("Could not find Module column. Using white unclassified strip.")
    print(colnames(obj))
    return(default_df)
  }
  
  class_col <- intersect(
    c(
      "ClassGroup",
      "class_group",
      "DominantClass",
      "dominant_class",
      "Dominant_Lipid_Class",
      "dominant_lipid_class",
      "LipidClass",
      "lipid_class",
      "Class",
      "dominant_hub_class",
      "DominantHubClass",
      "dominant_refmet_class"
    ),
    colnames(obj)
  )[1]
  
  if (is.na(class_col)) {
    message("Could not find lipid class column. Available columns:")
    print(colnames(obj))
    message("Using white unclassified strip.")
    return(default_df)
  }
  
  ann <- data.frame(
    Module     = normalize_module_names(obj$Module),
    LipidClass = as.character(obj[[class_col]]),
    stringsAsFactors = FALSE
  )
  
  ann <- ann[!duplicated(ann$Module), , drop = FALSE]
  
  ann$ClassGroup <- simplify_lipid_class(ann$LipidClass)
  ann$ClassColor <- unname(LIPID_CLASS_COLORS[ann$ClassGroup])
  ann$ClassColor[is.na(ann$ClassColor)] <- unname(LIPID_CLASS_COLORS["Unknown"])
  
  out <- merge(
    default_df[, "Module", drop = FALSE],
    ann,
    by = "Module",
    all.x = TRUE,
    sort = FALSE
  )
  
  out$LipidClass[is.na(out$LipidClass)] <- "Unknown"
  out$ClassGroup[is.na(out$ClassGroup)] <- "Unknown"
  out$ClassColor[is.na(out$ClassColor)] <- unname(LIPID_CLASS_COLORS["Unknown"])
  
  # Treat WGCNA grey module as its own category, not as a lipid class.
  grey_idx <- tolower(gsub("^ME", "", out$Module)) == "grey"
  out$LipidClass[grey_idx] <- "Grey module"
  out$ClassGroup[grey_idx] <- "Grey module"
  out$ClassColor[grey_idx] <- unname(LIPID_CLASS_COLORS["Grey module"])
  
  message("  Module lipid-class annotations loaded: ", sum(out$ClassGroup != "Unknown"))
  print(table(out$ClassGroup))
  
  out
}

make_final_module_order <- function(current_modules,
                                    module_order_file_order = NULL,
                                    module_class_df = NULL) {
  
  current_modules <- normalize_module_names(current_modules)
  
  if (is.null(module_order_file_order)) {
    module_order_file_order <- current_modules
  } else {
    module_order_file_order <- normalize_module_names(module_order_file_order)
  }
  
  module_order_file_order <- module_order_file_order[module_order_file_order %in% current_modules]
  remaining <- current_modules[!(current_modules %in% module_order_file_order)]
  base_order <- c(module_order_file_order, remaining)
  base_order <- base_order[!duplicated(base_order)]
  
  if (is.null(module_class_df)) {
    return(base_order)
  }
  
  class_map <- setNames(module_class_df$ClassGroup, module_class_df$Module)
  class_vec <- class_map[base_order]
  class_vec[is.na(class_vec)] <- "Unknown"
  
  class_rank <- match(class_vec, LIPID_CLASS_PRIORITY)
  class_rank[is.na(class_rank)] <- length(LIPID_CLASS_PRIORITY)
  
  final_order <- base_order[order(class_rank, seq_along(base_order))]
  
  message("Final module order places lipid classes in this order:")
  print(table(factor(class_vec[order(class_rank, seq_along(base_order))],
                     levels = LIPID_CLASS_PRIORITY)))
  
  final_order
}

# ── Load config and data ─────────────────────────────────

cfg <- yaml::read_yaml("config/config.yml")

method <- cfg$module_trait$method %||% "spearman"
method <- tolower(method)

wgcna_run_tag <- get_current_wgcna_run_tag(cfg)

message("WGCNA run: ", wgcna_run_tag, " | Method: ", toupper(method))

out <- file.path(cfg$output$s08, wgcna_run_tag, method)
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

# ── Module filter for Figure 2 ───────────────────────────
# Keeps only modules that have at least one significant
# association (p < cutoff OR q < cutoff) across ANY trait
# in ANY of the three timepoints.
# filter_mode: "nominal" | "fdr" | "none"

filter_modules_for_figure2 <- function(module_order,
                                       moduleTraitPvalue_ID,
                                       moduleTraitQvalue_ID,
                                       traits_use,
                                       filter_mode,
                                       cutoff_nominal,
                                       cutoff_fdr) {

  # ── Always remove the WGCNA grey module ──────────────────
  # Grey = unassigned lipids; not a meaningful biological module
  grey_idx <- tolower(gsub("^ME", "", module_order)) == "grey"
  if (any(grey_idx)) {
    message("  Removing WGCNA grey module from Figure 2 (unassigned lipids)")
  }
  module_order <- module_order[!grey_idx]

  if (filter_mode == "none") {
    message("  Module filter: none — showing all ", length(module_order), " modules (grey excluded)")
    return(module_order)
  }

  # Build a single logical matrix: rows = modules, cols = timepoints
  # TRUE means "this module has at least one significant trait at this timepoint"
  n_sets <- length(moduleTraitPvalue_ID)

  sig_any <- rep(FALSE, length(module_order))
  names(sig_any) <- module_order

  for (set in seq_len(n_sets)) {

    if (filter_mode == "nominal") {
      mat <- moduleTraitPvalue_ID[[set]]
      cutoff <- cutoff_nominal
    } else if (filter_mode == "fdr") {
      mat <- moduleTraitQvalue_ID[[set]]
      cutoff <- cutoff_fdr
    } else {
      stop("Unknown figure2_module_filter value: '", filter_mode,
           "'. Use 'nominal', 'fdr', or 'none'.")
    }

    # Restrict to the traits shown in Figure 2
    traits_in_mat <- intersect(traits_use, colnames(mat))
    if (length(traits_in_mat) == 0) next

    mat_sub <- mat[, traits_in_mat, drop = FALSE]

    # Normalise row names so they match module_order (ME-prefixed)
    rownames(mat_sub) <- normalize_module_names(rownames(mat_sub))

    modules_in_mat <- intersect(module_order, rownames(mat_sub))
    if (length(modules_in_mat) == 0) next

    # Any trait below cutoff for this timepoint?
    hits <- apply(mat_sub[modules_in_mat, , drop = FALSE], 1,
                  function(x) any(!is.na(x) & x < cutoff))

    sig_any[modules_in_mat] <- sig_any[modules_in_mat] | hits
  }

  filtered <- module_order[sig_any]

  message(sprintf(
    "  Module filter (%s, cutoff = %s): %d / %d modules kept for Figure 2",
    filter_mode,
    ifelse(filter_mode == "fdr", cutoff_fdr, cutoff_nominal),
    length(filtered),
    length(module_order)
  ))
  message("  Removed modules (go to supplementary): ",
          paste(gsub("^ME", "", setdiff(module_order, filtered)), collapse = ", "))

  filtered
}

# ── Full all-trait heatmap function ──────────────────────
# Keeps WGCNA::labeledHeatmap for supplementary/full heatmaps.
# No lipid-class strip here because labeledHeatmap cannot safely show
# both row color annotation and module labels.

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
  
  if (!is.null(module_order)) {
    keep_order <- module_order[module_order %in% rownames(cor_mat)]
    remaining  <- rownames(cor_mat)[!(rownames(cor_mat) %in% keep_order)]
    final_order <- c(keep_order, remaining)
    
    cor_mat  <- cor_mat[final_order, , drop = FALSE]
    pval_mat <- pval_mat[final_order, , drop = FALSE]
    qval_mat <- qval_mat[final_order, , drop = FALSE]
  }
  
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
  
  panel_width_mm  <- FIG_WIDTH_FULL
  panel_height_mm <- max(175, nrow(cor_mat) * 5.8)
  panel_height_mm <- min(panel_height_mm, 220)
  
  heatmap_cex_text <- 0.55
  heatmap_cex_y    <- 0.62
  heatmap_cex_x    <- 0.68
  
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

# ── Direct journal-quality Figure 2: ggplot version ──────
# This version supports:
#   1. True block-style lipid-class strip
#   2. Preserved module labels
#   3. Separate lipid-class legend
#   4. Configurable panel labels: lower-case or upper-case
#   5. Shared x-axis labels shown only on the bottom panel
#   6. BH/FDR note in the legend area, not repeated in panels

make_direct_figure2_heatmaps <- function(moduleTraitCor_ID,
                                         moduleTraitPvalue_ID,
                                         moduleTraitQvalue_ID,
                                         module_order = NULL,
                                         module_class_df = NULL,
                                         cfg,
                                         out_dir) {
  
  message("\nCreating direct journal-quality Figure 2 heatmaps with lipid-class block strip...")
  
  trait_set    <- cfg$module_trait$figure2_trait_set %||% "primary"
  heatmap_mode <- cfg$module_trait$figure2_heatmap_mode %||% "fdr"
  
  panel_label_case <- tolower(cfg$module_trait$figure2_panel_label_case %||% "lower")
  shared_x_axis    <- cfg$module_trait$figure2_shared_x_axis %||% TRUE
  fig_height_mm    <- cfg$module_trait$figure2_height_mm %||% FIG_HEIGHT_MAX
  
  if (!(panel_label_case %in% c("lower", "upper"))) {
    warning("Invalid figure2_panel_label_case. Using 'lower'.")
    panel_label_case <- "lower"
  }
  
  panel_letters <- if (panel_label_case == "upper") {
    c("A", "B", "C")
  } else {
    c("a", "b", "c")
  }
  
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
  message("  Panel label case: ", panel_label_case)
  message("  Shared x-axis: ", shared_x_axis)
  message("  Figure height: ", fig_height_mm, " mm")
  
  fig_dir <- file.path("results", "figures", "final")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  
  fig_pdf <- file.path(fig_dir, paste0("Figure2_direct_", trait_set, "_", heatmap_mode, ".pdf"))
  fig_png <- file.path(fig_dir, paste0("Figure2_direct_", trait_set, "_", heatmap_mode, ".png"))
  fig_tif <- file.path(fig_dir, paste0("Figure2_direct_", trait_set, "_", heatmap_mode, ".tiff"))
  
  panel_titles <- c("10–16 weeks", "36–38 weeks", "Postpartum")
  
  if (is.null(module_order)) {
    module_order <- rownames(moduleTraitCor_ID[[1]])
  }
  
  module_order <- normalize_module_names(module_order)
  
  if (is.null(module_class_df)) {
    module_class_df <- data.frame(
      Module = module_order,
      ClassGroup = "Unknown",
      ClassColor = unname(LIPID_CLASS_COLORS["Unknown"]),
      stringsAsFactors = FALSE
    )
  }
  
  module_class_df$Module <- normalize_module_names(module_class_df$Module)
  
  # Ensure grey module is its own category.
  grey_idx <- tolower(gsub("^ME", "", module_class_df$Module)) == "grey"
  module_class_df$ClassGroup[grey_idx] <- "Grey module"
  module_class_df$ClassColor[grey_idx] <- unname(LIPID_CLASS_COLORS["Grey module"])
  
  current_modules <- normalize_module_names(rownames(moduleTraitCor_ID[[1]]))
  module_order <- module_order[module_order %in% current_modules]
  missing_modules <- current_modules[!(current_modules %in% module_order)]
  module_order <- c(module_order, missing_modules)
  module_order <- module_order[!duplicated(module_order)]
  
  # ── Filter modules for Figure 2 ────────────────────────
  filter_mode      <- tolower(cfg$module_trait$figure2_module_filter %||% "none")
  cutoff_nominal   <- cfg$module_trait$figure2_module_filter_cutoff_nominal %||% raw_p_cutoff
  cutoff_fdr_fig   <- cfg$module_trait$figure2_module_filter_cutoff_fdr     %||% fdr_cutoff

  module_order <- filter_modules_for_figure2(
    module_order            = module_order,
    moduleTraitPvalue_ID    = moduleTraitPvalue_ID,
    moduleTraitQvalue_ID    = moduleTraitQvalue_ID,
    traits_use              = traits_use,
    filter_mode             = filter_mode,
    cutoff_nominal          = cutoff_nominal,
    cutoff_fdr              = cutoff_fdr_fig
  )

  if (length(module_order) == 0) {
    warning("No modules passed the Figure 2 filter. Check filter settings in config.yml.")
    return(invisible(NULL))
  }
  # ── End module filter ───────────────────────────────────

  make_panel_plot <- function(set_index) {
    
    show_x_axis <- if (isTRUE(shared_x_axis)) {
      set_index == length(moduleTraitCor_ID)
    } else {
      TRUE
    }
    
    cor_mat <- moduleTraitCor_ID[[set_index]]
    p_mat   <- moduleTraitPvalue_ID[[set_index]]
    q_mat   <- moduleTraitQvalue_ID[[set_index]]
    
    cor_mat <- cor_mat[, traits_use, drop = FALSE]
    p_mat   <- p_mat[, traits_use, drop = FALSE]
    q_mat   <- q_mat[, traits_use, drop = FALSE]
    
    rownames(cor_mat) <- normalize_module_names(rownames(cor_mat))
    rownames(p_mat)   <- normalize_module_names(rownames(p_mat))
    rownames(q_mat)   <- normalize_module_names(rownames(q_mat))
    
    keep_order <- module_order[module_order %in% rownames(cor_mat)]
    
    cor_mat <- cor_mat[keep_order, , drop = FALSE]
    p_mat   <- p_mat[keep_order, , drop = FALSE]
    q_mat   <- q_mat[keep_order, , drop = FALSE]
    
    if (heatmap_mode == "fdr") {
      sig_mat <- q_mat
      cutoff <- fdr_cutoff
    } else {
      sig_mat <- p_mat
      cutoff <- raw_p_cutoff
    }
    
    # Numeric y positions. First module in keep_order appears at top.
    module_y <- setNames(rev(seq_along(keep_order)), keep_order)
    
    # Heatmap data.
    df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
    colnames(df) <- c("Module", "Trait", "Correlation")
    
    sig_df <- as.data.frame(as.table(sig_mat), stringsAsFactors = FALSE)
    colnames(sig_df) <- c("Module", "Trait", "SigValue")
    
    df <- dplyr::left_join(df, sig_df, by = c("Module", "Trait"))
    
    df$Module <- normalize_module_names(df$Module)
    df$Trait  <- as.character(df$Trait)
    df$Star <- ifelse(!is.na(df$SigValue) & df$SigValue < cutoff, "*", "")
    
    trait_positions <- seq_along(traits_use)
    names(trait_positions) <- traits_use
    
    df$TraitX  <- trait_positions[df$Trait]
    df$ModuleY <- module_y[df$Module]
    
    # Annotation data matched to ordered modules.
    ann_df <- module_class_df[, c("Module", "ClassGroup", "ClassColor"), drop = FALSE]
    ann_df$Module <- normalize_module_names(ann_df$Module)
    
    ann_ordered <- data.frame(Module = keep_order, stringsAsFactors = FALSE)
    ann_ordered <- dplyr::left_join(ann_ordered, ann_df, by = "Module")
    
    ann_ordered$ClassGroup[is.na(ann_ordered$ClassGroup)] <- "Unknown"
    ann_ordered$ClassColor[is.na(ann_ordered$ClassColor)] <- unname(LIPID_CLASS_COLORS["Unknown"])
    
    grey_idx2 <- tolower(gsub("^ME", "", ann_ordered$Module)) == "grey"
    ann_ordered$ClassGroup[grey_idx2] <- "Grey module"
    ann_ordered$ClassColor[grey_idx2] <- unname(LIPID_CLASS_COLORS["Grey module"])
    
    ann_ordered$ModuleY <- module_y[ann_ordered$Module]
    
    # Create block-level class segments for contiguous modules.
    r <- rle(ann_ordered$ClassGroup)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    
    class_blocks <- data.frame(
      ClassGroup = r$values,
      start = starts,
      end = ends,
      stringsAsFactors = FALSE
    )
    
    class_blocks$ClassColor <- ann_ordered$ClassColor[class_blocks$start]
    
    class_blocks$ymin <- pmin(
      ann_ordered$ModuleY[class_blocks$start],
      ann_ordered$ModuleY[class_blocks$end]
    ) - 0.5
    
    class_blocks$ymax <- pmax(
      ann_ordered$ModuleY[class_blocks$start],
      ann_ordered$ModuleY[class_blocks$end]
    ) + 0.5
    
    # Abbreviated names for strip labels
    abbrev_map <- c(
        "TG"          = "TG",
        "DG"          = "DG",
        "CE"          = "CE",
        "PC/LPC"      = "PC/LPC",
        "PE/LPE"      = "PE/LPE",
        "SM"          = "SM",
        "Cer"         = "Cer",
        "FA"          = "FA",
        "Oth"       = "Oth",
        "Unknown"     = "",
        "Grey module" = ""        # grey excluded already but safe fallback
      )

      class_blocks$ClassAbbrev <- unname(abbrev_map[class_blocks$ClassGroup])
      class_blocks$ClassAbbrev[is.na(class_blocks$ClassAbbrev)] <- ""

      # White text on dark backgrounds, black text on light backgrounds
      # Dark classes: TG (red), PC/LPC (dark blue), SM (teal), FA (brown)
      dark_classes <- c("TG", "PC/LPC", "SM", "FA")
      class_blocks$TextColor <- ifelse(
        class_blocks$ClassGroup %in% dark_classes,
        "white",
        "black"
      )

    p <- ggplot2::ggplot() +
      
      # Block-style lipid-class strip.
      ggplot2::geom_rect(
        data = class_blocks,
        ggplot2::aes(
          xmin = 0.10,
          xmax = 0.52, # widened from 0.34
          ymin = ymin,
          ymax = ymax
        ),
        fill = class_blocks$ClassColor,
        colour = "black",
        linewidth = 0.12,
        inherit.aes = FALSE
      ) +

      # Abbreviated lipid class labels inside the strip
      ggplot2::geom_text(
        data = class_blocks,
        ggplot2::aes(
          x     = 0.31,           # horizontal centre of the strip
          y     = (ymin + ymax) / 2,
          label = ClassAbbrev,
          colour = TextColor
        ),
        angle    = 90,
        size     = 2.3,           # ~6.5pt — fits inside narrow blocks
        fontface = "bold",
        family   = "Helvetica",
        inherit.aes = FALSE
      ) +
      ggplot2::scale_colour_identity() +   # tells ggplot colour values are literal
      
        
      # Correlation heatmap.
      ggplot2::geom_tile(
        data = df,
        ggplot2::aes(x = TraitX, y = ModuleY, fill = Correlation),
        colour = "white",
        linewidth = 0.12,
        width = 0.98,
        height = 0.95
      ) +
      
      ggplot2::annotate(
        "text",
        x      = -6.3,                         # sits left of the heatmap # was -1.2, now further left
        y      = mean(range(df$ModuleY)),       # vertically centred on heatmap
        label  = paste0(panel_titles[set_index]),
        angle  = 90,
        hjust  = 0.5,
        vjust  = 1,
        size   = 5.2,
        colour   = "black",                          # ~9pt
        # fontface = "bold",
        family = "Helvetica"
      ) +

      # Significance stars.
      ggplot2::geom_text(
        data = df[df$Star != "", , drop = FALSE],
        ggplot2::aes(x = TraitX, y = ModuleY, label = Star),
        size = 2.8,
        colour = "black",
        hjust   = 0.5,        # horizontal centre
        vjust   = 0.5, 
        family = "Helvetica"
      ) +
      
      # WGCNA-style blue-white-red heatmap scale.
      ggplot2::scale_fill_gradientn(
        colors = WGCNA::blueWhiteRed(100),
        limits = c(-1, 1),
        breaks = c(-1, -0.5, 0, 0.5, 1),
        name = ifelse(method == "pearson", "Pearson r", "Spearman \u03c1")
      ) +
      
      ggplot2::scale_x_continuous(
        breaks = trait_positions,
        labels = if (show_x_axis) traits_use else rep("", length(traits_use)),
        expand = ggplot2::expansion(mult = c(0.01, 0.02))
      ) +
      
      ggplot2::scale_y_continuous(
        breaks = module_y[keep_order],
        labels = gsub("^ME", "", keep_order),
        expand = ggplot2::expansion(mult = c(0.002, 0.002))
      ) +
      
      ggplot2::coord_cartesian(
        xlim = c(0.05, length(traits_use) + 0.5),
        clip = "off"
      ) +
      
      ggplot2::labs(
        title = NULL,                          # title moved to left margin
        tag   = panel_letters[set_index],
        x     = NULL,
        y     = NULL
        ) +

      ggplot2::theme_minimal(base_family = "Helvetica", base_size = 8) +
       # Panel title as rotated text in left margin
  

      ggplot2::theme(
        plot.title = ggplot2::element_blank(),  # suppressed — using annotate instead
        plot.tag   = ggplot2::element_text(
          face = "bold",
          size = 13
        ),
        plot.tag.position = c(-0.085, 1.0),    # slightly lower since no title above
        axis.text.x = ggplot2::element_text(
          angle    = 45, # was 45
          hjust    = 1,
          vjust    = 1, # change from 1 to 0.5 — centers labels vertically at tick
          size     = if (show_x_axis) 8 else 0,
          colour   = "black",
          face     = "bold"
        ),
        axis.ticks.x = if (show_x_axis) {
          ggplot2::element_line(linewidth = 0.25)
        } else {
          ggplot2::element_blank()
        },
        axis.text.y = ggplot2::element_text(
          size   = 7, #was 8, reduced to fit longer module labels with class strip
          colour = "black",
          face   = "bold"
        ),
        panel.grid   = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(
          colour    = "black",
          fill      = NA,
          linewidth = 0.4
        ),
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 8, face = "bold"),
        legend.text  = ggplot2::element_text(size = 8, face = "bold"),
        plot.margin  = ggplot2::margin(
          t = 4,                               # reduced top — no title above
          r = 14,
          b = if (show_x_axis) 8 else 1,
          l = 95                               # increased left: module labels + rotated title # was 72, increased to accommodate title position
        )
      )
      
    p
  }
  
  make_class_legend_plot <- function(module_class_df) {
    
    present <- unique(module_class_df$ClassGroup)
    present <- present[!is.na(present)]
    present <- present[present %in% names(LIPID_CLASS_COLORS)]
    present <- LIPID_CLASS_PRIORITY[LIPID_CLASS_PRIORITY %in% present]
    
    if (length(present) == 0) {
      present <- "Unknown"
    }
    
    legend_labels <- c(
      "TG"          = "Triglycerides (TG)",
      "DG"          = "Diacylglycerols (DG)",
      "CE"          = "Cholesteryl esters (CE)",
      "PC/LPC"      = "PC/LPC",
      "PE/LPE"      = "PE/LPE",
      "SM"          = "Sphingomyelins (SM)",
      "Cer"         = "Ceramides (Cer)",
      "FA"          = "Fatty acyls (FA)",
      "Oth"       = "Other class",
      "Unknown"     = "Unclassified",
      "Grey module" = "WGCNA grey module"
    )
    
    leg_df <- data.frame(
      ClassGroup = present,
      Label = unname(legend_labels[present]),
      Color = unname(LIPID_CLASS_COLORS[present]),
      stringsAsFactors = FALSE
    )
    
    n_col <- 4
    leg_df$col <- ((seq_len(nrow(leg_df)) - 1) %% n_col) + 1
    leg_df$row <- ceiling(seq_len(nrow(leg_df)) / n_col)
    max_row <- max(leg_df$row)
    leg_df$y <- max_row - leg_df$row + 1
    
    sig_note <- if (heatmap_mode == "fdr") {
      paste0("* family-based BH q < ", fdr_cutoff)
    } else {
      paste0("* nominal p < ", raw_p_cutoff)
    }
    
    sig_note_df <- data.frame(
      x = 0.8,
      y = -0.1,
      label = sig_note,
      stringsAsFactors = FALSE
    )
    
    ggplot2::ggplot(leg_df) +
      ggplot2::geom_tile(
        ggplot2::aes(x = col, y = y),
        fill = leg_df$Color,
        colour = "black",
        width = 0.12,
        height = 0.35
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = col + 0.10, y = y, label = Label),
        hjust = 0,
        size = 2.8,               # slight bump to ~8pt equivalent
        family = "Helvetica",
        fontface = "bold"         # ADD
      ) +
      ggplot2::geom_text(
        data = sig_note_df,
        ggplot2::aes(x = x, y = y, label = label),
        hjust = 0,
        size = 2.8,
        family = "Helvetica",
        fontface = "bold",        # ADD
        inherit.aes = FALSE
      ) +
      ggplot2::xlim(0.7, n_col + 1.7) +
      ggplot2::ylim(-0.3, max_row + 0.6) +
      ggplot2::theme_void(base_family = "Helvetica") +
      ggplot2::theme(
        plot.margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 40)
      )
  }
  
  p1 <- make_panel_plot(1)
  p2 <- make_panel_plot(2)
  p3 <- make_panel_plot(3)
  # p_legend <- make_class_legend_plot(module_class_df)
  
  message(
    "  Class strip: block colors represent dominant lipid class; ",
    "grey is reserved only for the WGCNA grey module."
  )
  
  draw_combined <- function() {
    grid::grid.newpage()
    
    grid::pushViewport(
      grid::viewport(
        layout = grid::grid.layout(
          nrow = 3, #was 4, now 3 since legend is separate
          ncol = 1,
          heights = grid::unit(c(1.0, 1.0, 1.1), "null")
        )
      )
    )
    
    print(p1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p2, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(p3, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
    
    grid::popViewport()
  }
  
  fig_width_mm <- FIG_WIDTH_FULL
  
  convert_pdf_to_rasters <- function(pdf_file, png_file, tiff_file, dpi = 600) {
    
    message("Converting PDF to high-resolution PNG/TIFF...")
    
    pdftocairo_path <- Sys.which("pdftocairo")
    
    if (pdftocairo_path == "") {
      warning(
        "pdftocairo not found. Skipping PNG/TIFF conversion.\n",
        "Install poppler or use the PDF directly."
      )
      return(invisible(NULL))
    }
    
    tmp_prefix <- tempfile("figure2_pdf_render_")
    
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
  
  grDevices::cairo_pdf(
    filename = fig_pdf,
    width = fig_width_mm / 25.4,
    height = fig_height_mm / 25.4,
    family = "Helvetica"
  )
  draw_combined()
  dev.off()
  
  convert_pdf_to_rasters(
    pdf_file  = fig_pdf,
    png_file  = fig_png,
    tiff_file = fig_tif,
    dpi       = FIGURE_DPI_FINAL
  )
  
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

# Load module order and lipid-class annotation.
first_me <- consensusMEs_ID[[1]]$data

module_order_file_order <- load_module_order_file(
  cfg = cfg,
  wgcna_run_tag = wgcna_run_tag,
  current_modules = colnames(first_me)
)

module_class_df <- load_module_class_annotation(
  cfg = cfg,
  wgcna_run_tag = wgcna_run_tag,
  current_modules = colnames(first_me)
)

module_order <- make_final_module_order(
  current_modules = colnames(first_me),
  module_order_file_order = module_order_file_order,
  module_class_df = module_class_df
)

# Save exact module class/order used.
module_class_order_used <- module_class_df[match(module_order, module_class_df$Module), , drop = FALSE]

write.csv(
  module_class_order_used,
  file.path(out, "module_class_order_used_in_heatmaps.csv"),
  row.names = FALSE
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
  
  message(
    "  FDR q<", fdr_cutoff, ": ",
    sum(q_mat < fdr_cutoff, na.rm = TRUE),
    " module-trait cells"
  )
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

# Save Excel files with correlations, p-values, q-values.
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
  module_class_df      = module_class_df,
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

message("\nScript 08 complete → ", out)
message("Heatmap folders created:")
for (hf in unique(heatmap_folders)) {
  message("  - ", hf)
}
message("Figure 2 panels saved to results/figures/panels/fig2/<heatmap_mode>/")
message("Direct Figure 2 saved to results/figures/final/")
message("Module class/order file saved to: ", file.path(out, "module_class_order_used_in_heatmaps.csv"))

# # R/08_module_trait_correlations.R
# # -------------------------------------------------------
# # Module-trait correlations with configurable raw-p/FDR heatmaps
# # Method: pearson/spearman from config.yml
# #
# # Outputs:
# #   results/05_module_trait_cor/<run_tag>/<method>/
# #   results/figures/panels/fig2/<heatmap_mode>/
# #   results/figures/final/Figure2_direct_<trait_set>_<mode>.pdf/png/tiff
# # -------------------------------------------------------

# source("R/00_load_packages.R")
# source("R/00_figure_theme.R")

# # ── Helpers ──────────────────────────────────────────────
# `%||%` <- function(a, b) {
#   if (!is.null(a)) a else b
# }

# get_current_wgcna_run_tag <- function(cfg) {
#   tag_file <- file.path(cfg$output$processed, "current_wgcna_run_tag.txt")
  
#   if (length(tag_file) != 1 || is.na(tag_file) || tag_file == "") {
#     stop("Invalid current_wgcna_run_tag.txt path. Check cfg$output$processed.")
#   }
  
#   if (!file.exists(tag_file)) {
#     stop(
#       "Missing current WGCNA run tag file: ", tag_file, "\n",
#       "Run R/03_run_consensus_wgcna.R first."
#     )
#   }
  
#   trimws(readLines(tag_file, warn = FALSE)[1])
# }

# cutoff_label <- function(x) {
#   # Gives common folder names:
#   # 0.05 -> 005
#   # 0.10 -> 01
#   if (abs(x - 0.05) < 1e-9) return("005")
#   if (abs(x - 0.10) < 1e-9) return("01")
#   out <- sprintf("%.3f", x)
#   out <- sub("^0\\.", "", out)
#   out <- sub("0+$", "", out)
#   out
# }

# normalize_module_names <- function(x) {
#   x <- as.character(x)
#   x <- gsub("^ME", "", x)
#   paste0("ME", x)
# }

# load_module_order <- function(cfg, wgcna_run_tag, current_modules) {
  
#   module_order_file <- file.path(
#     cfg$output$processed,
#     paste0("module_order_by_hub_lipid_class_", wgcna_run_tag, ".rds")
#   )
  
#   if (!file.exists(module_order_file)) {
#     message("No hub-lipid-class module order file found. Using current WGCNA module order.")
#     return(current_modules)
#   }
  
#   message("Loading module order from: ", module_order_file)
#   obj <- readRDS(module_order_file)
  
#   if (is.data.frame(obj)) {
#     if ("Module" %in% colnames(obj)) {
#       module_order <- obj$Module
#     } else if ("module" %in% colnames(obj)) {
#       module_order <- obj$module
#     } else {
#       message("Module order file has no Module/module column. Using current order.")
#       return(current_modules)
#     }
#   } else {
#     module_order <- obj
#   }
  
#   module_order <- normalize_module_names(module_order)
#   current_modules_norm <- normalize_module_names(current_modules)
  
#   ordered_keep <- module_order[module_order %in% current_modules_norm]
#   remaining <- current_modules_norm[!(current_modules_norm %in% ordered_keep)]
  
#   final_order <- c(ordered_keep, remaining)
#   final_order <- final_order[!duplicated(final_order)]
  
#   message("  Ordered modules loaded: ", length(ordered_keep))
#   message("  Remaining modules appended: ", length(remaining))
  
#   final_order
# }

# # ── Load config and data ─────────────────────────────────
# cfg <- yaml::read_yaml("config/config.yml")

# method <- cfg$module_trait$method %||% "spearman"
# method <- tolower(method)

# wgcna_run_tag <- get_current_wgcna_run_tag(cfg)

# message("WGCNA run: ", wgcna_run_tag, " | Method: ", toupper(method))

# out <- file.path(cfg$output$s08, wgcna_run_tag, method)
# dir.create(out, showWarnings = FALSE, recursive = TRUE)

# consensusMEs_ID <- readRDS(file.path(
#   cfg$output$processed,
#   paste0("consensusMEs_ID_", wgcna_run_tag, ".rds")
# ))

# demographic_data <- readRDS(file.path(
#   cfg$output$processed,
#   "demographic_data_bmi.rds"
# ))

# exprSize_ID <- readRDS(file.path(
#   cfg$output$processed,
#   "exprSize_ID.rds"
# ))

# All_nSets_ID <- exprSize_ID$nSets

# tp_labels  <- c("Baseline", "Wk36_38", "Postpartum")
# tp_display <- c("10–16 weeks", "36–38 weeks", "Postpartum")
# panel_ids  <- c("A", "B", "C")

# # Keep only numeric phenotype traits for correlation
# pheno_data <- demographic_data
# pheno_data <- pheno_data[, sapply(pheno_data, is.numeric), drop = FALSE]

# # Heatmap / FDR options
# heatmap_modes <- unlist(cfg$module_trait$heatmap_significance %||% "raw")
# heatmap_modes <- tolower(heatmap_modes)

# raw_p_cutoff <- cfg$module_trait$raw_p_cutoff %||% 0.05
# fdr_cutoff   <- cfg$module_trait$fdr_cutoff %||% 0.05
# fdr_method   <- cfg$module_trait$fdr_method %||% "BH"
# fdr_scope    <- cfg$module_trait$fdr_scope %||% "family"

# message("Heatmap significance modes: ", paste(heatmap_modes, collapse = ", "))
# message(
#   "Raw p cutoff: ", raw_p_cutoff,
#   " | FDR cutoff: ", fdr_cutoff,
#   " | FDR method: ", fdr_method,
#   " | FDR scope: ", fdr_scope
# )

# # ── Correlation function ─────────────────────────────────
# compute_cor <- function(x, y, method) {
  
#   if (method == "spearman") {
    
#     cor_mat <- stats::cor(
#       x,
#       y,
#       method = "spearman",
#       use = "pairwise.complete.obs"
#     )
    
#     n <- nrow(x)
    
#     # Approximate p-values using t approximation
#     # Protect against r = +/- 1
#     cor_for_test <- pmin(pmax(cor_mat, -0.999999), 0.999999)
#     t_stat <- cor_for_test * sqrt((n - 2) / (1 - cor_for_test^2))
#     pval_mat <- 2 * stats::pt(-abs(t_stat), df = n - 2)
    
#   } else if (method == "pearson") {
    
#     cor_mat <- WGCNA::cor(
#       x,
#       y,
#       method = "pearson",
#       use = "p"
#     )
    
#     pval_mat <- WGCNA::corPvalueFisher(
#       cor_mat,
#       nrow(x),
#       twoSided = TRUE
#     )
    
#   } else {
#     stop("Unsupported correlation method: ", method)
#   }
  
#   dimnames(pval_mat) <- dimnames(cor_mat)
  
#   list(cor = cor_mat, pval = pval_mat)
# }

# # ── Trait families for FDR ───────────────────────────────
# make_trait_family_vector <- function(traits, cfg) {
  
#   families <- cfg$module_trait$fdr_families
  
#   trait_family <- rep("Exploratory_other", length(traits))
#   names(trait_family) <- traits
  
#   if (!is.null(families)) {
#     for (fam in names(families)) {
#       fam_traits <- unlist(families[[fam]])
#       fam_traits <- fam_traits[fam_traits %in% traits]
#       trait_family[fam_traits] <- fam
#     }
#   }
  
#   trait_family
# }

# apply_fdr <- function(p_mat, trait_family, method = "BH", scope = "family") {
  
#   q_mat <- matrix(
#     NA_real_,
#     nrow = nrow(p_mat),
#     ncol = ncol(p_mat),
#     dimnames = dimnames(p_mat)
#   )
  
#   if (scope == "all_traits") {
    
#     q_mat[] <- p.adjust(as.vector(p_mat), method = method)
    
#   } else if (scope == "family") {
    
#     for (fam in unique(trait_family)) {
      
#       fam_traits <- names(trait_family)[trait_family == fam]
#       fam_traits <- fam_traits[fam_traits %in% colnames(p_mat)]
      
#       if (length(fam_traits) == 0) next
      
#       p_sub <- p_mat[, fam_traits, drop = FALSE]
      
#       q_sub <- matrix(
#         p.adjust(as.vector(p_sub), method = method),
#         nrow = nrow(p_sub),
#         ncol = ncol(p_sub),
#         dimnames = dimnames(p_sub)
#       )
      
#       q_mat[, fam_traits] <- q_sub
#     }
    
#   } else {
#     stop("Unsupported fdr_scope: ", scope)
#   }
  
#   q_mat
# }

# # ── Heatmap function for all-trait panels ─────────────────
# make_heatmap_panel <- function(cor_mat,
#                                pval_mat,
#                                qval_mat,
#                                title_str,
#                                pdf_path,
#                                fig_panel_id,
#                                heatmap_mode = c("raw", "fdr"),
#                                module_order = NULL) {
  
#   heatmap_mode <- match.arg(heatmap_mode)
  
#   if (is.null(cor_mat) || nrow(cor_mat) == 0 || ncol(cor_mat) == 0) {
#     return(invisible(NULL))
#   }
  
#   # Reorder modules, with TG-associated modules first if module_order exists
#   if (!is.null(module_order)) {
#     keep_order <- module_order[module_order %in% rownames(cor_mat)]
#     remaining  <- rownames(cor_mat)[!(rownames(cor_mat) %in% keep_order)]
#     final_order <- c(keep_order, remaining)
    
#     cor_mat  <- cor_mat[final_order, , drop = FALSE]
#     pval_mat <- pval_mat[final_order, , drop = FALSE]
#     qval_mat <- qval_mat[final_order, , drop = FALSE]
#   }
  
#   # Significance matrix
#   if (heatmap_mode == "fdr") {
#     sig_mat <- qval_mat
#     cutoff  <- fdr_cutoff
#     sig_note <- paste0("BH q < ", cutoff)
#   } else {
#     sig_mat <- pval_mat
#     cutoff  <- raw_p_cutoff
#     sig_note <- paste0("p < ", cutoff)
#   }
  
#   textMatrix <- apply(sig_mat, 2, function(x) {
#     sapply(x, function(v) {
#       ifelse(!is.na(v) && v < cutoff, "*", "")
#     })
#   })
  
#   dim(textMatrix) <- dim(cor_mat)
  
#   row_labels <- gsub("^ME", "", rownames(cor_mat))
  
#   # Dynamic panel size
#   panel_width_mm  <- FIG_WIDTH_FULL
#   panel_height_mm <- max(175, nrow(cor_mat) * 5.8)
#   panel_height_mm <- min(panel_height_mm, 220)
  
#   heatmap_cex_text <- 0.55
#   heatmap_cex_y    <- 0.62
#   heatmap_cex_x    <- 0.68
  
#   # Save vector PDF
#   grDevices::cairo_pdf(
#     filename = pdf_path,
#     width    = panel_width_mm / 25.4,
#     height   = panel_height_mm / 25.4,
#     family   = "Helvetica"
#   )
  
#   par(
#     mar = c(8.5, 7.8, 1.8, 3.5),
#     oma = c(0, 0, 0, 0),
#     family = "Helvetica",
#     xpd = NA
#   )
  
#   WGCNA::labeledHeatmap(
#     Matrix        = cor_mat,
#     xLabels       = colnames(cor_mat),
#     yLabels       = rownames(cor_mat),
#     ySymbols      = row_labels,
#     colorLabels   = FALSE,
#     colors        = WGCNA::blueWhiteRed(50),
#     textMatrix    = textMatrix,
#     setStdMargins = FALSE,
#     cex.text      = heatmap_cex_text,
#     textAdj       = c(0.5, 0.8),
#     zlim          = c(-1, 1),
#     main          = title_str,
#     cex.lab.y     = heatmap_cex_y,
#     cex.lab.x     = heatmap_cex_x,
#     plotLegend    = TRUE,
#     legendLabel   = ifelse(method == "pearson", "Pearson r", "Spearman \u03c1")
#   )
  
#   mtext(sig_note, side = 1, line = 7.2, adj = 1, cex = 0.65)
  
#   dev.off()
  
#   # Save 600 dpi PNG
#   png_path <- gsub("\\.pdf$", ".png", pdf_path)
  
#   png(
#     filename = png_path,
#     width    = panel_width_mm,
#     height   = panel_height_mm,
#     units    = "mm",
#     res      = FIGURE_DPI_FINAL,
#     type     = "cairo",
#     family   = "Helvetica"
#   )
  
#   par(
#     mar = c(8.5, 7.8, 1.8, 3.5),
#     oma = c(0, 0, 0, 0),
#     family = "Helvetica",
#     xpd = NA
#   )
  
#   WGCNA::labeledHeatmap(
#     Matrix        = cor_mat,
#     xLabels       = colnames(cor_mat),
#     yLabels       = rownames(cor_mat),
#     ySymbols      = row_labels,
#     colorLabels   = FALSE,
#     colors        = WGCNA::blueWhiteRed(50),
#     textMatrix    = textMatrix,
#     setStdMargins = FALSE,
#     cex.text      = heatmap_cex_text,
#     textAdj       = c(0.5, 0.8),
#     zlim          = c(-1, 1),
#     main          = title_str,
#     cex.lab.y     = heatmap_cex_y,
#     cex.lab.x     = heatmap_cex_x,
#     plotLegend    = TRUE,
#     legendLabel   = ifelse(method == "pearson", "Pearson r", "Spearman \u03c1")
#   )
  
#   mtext(sig_note, side = 1, line = 7.2, adj = 1, cex = 0.65)
  
#   dev.off()
  
#   # Register panel in heatmap-specific panel folder
#   heatmap_folder_name <- basename(dirname(pdf_path))
#   panel_dir <- file.path("results", "figures", "panels", "fig2", heatmap_folder_name)
#   dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)
  
#   file.copy(
#     png_path,
#     file.path(panel_dir, paste0("fig2_", fig_panel_id, ".png")),
#     overwrite = TRUE
#   )
  
#   file.copy(
#     pdf_path,
#     file.path(panel_dir, paste0("fig2_", fig_panel_id, ".pdf")),
#     overwrite = TRUE
#   )
  
#   message(
#     "  ", heatmap_mode, " panel fig2_", fig_panel_id,
#     " saved | height: ", round(panel_height_mm, 1), " mm → ",
#     dirname(pdf_path)
#   )
  
#   invisible(list(pdf = pdf_path, png = png_path))
# }

# # ── Direct journal-quality Figure 2: primary traits only ──
# make_direct_figure2_heatmaps <- function(moduleTraitCor_ID,
#                                          moduleTraitPvalue_ID,
#                                          moduleTraitQvalue_ID,
#                                          module_order = NULL,
#                                          cfg,
#                                          out_dir) {
  
#   message("\nCreating direct journal-quality Figure 2 heatmaps...")
  
#   trait_set    <- cfg$module_trait$figure2_trait_set %||% "primary"
#   heatmap_mode <- cfg$module_trait$figure2_heatmap_mode %||% "fdr"
  
#   if (trait_set == "primary") {
#     traits_use <- unlist(cfg$module_trait$primary_figure2_traits)
#   } else {
#     traits_use <- colnames(moduleTraitCor_ID[[1]])
#   }
  
#   traits_use <- traits_use[traits_use %in% colnames(moduleTraitCor_ID[[1]])]
  
#   if (length(traits_use) == 0) {
#     warning("No direct Figure 2 traits found. Skipping direct Figure 2.")
#     return(invisible(NULL))
#   }
  
#   message("  Figure 2 trait set: ", trait_set)
#   message("  Traits used: ", paste(traits_use, collapse = ", "))
#   message("  Significance mode: ", heatmap_mode)
  
#   fig_dir <- file.path("results", "figures", "final")
#   dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  
#   fig_pdf <- file.path(fig_dir, paste0("Figure2_direct_", trait_set, "_", heatmap_mode, ".pdf"))
#   fig_png <- file.path(fig_dir, paste0("Figure2_direct_", trait_set, "_", heatmap_mode, ".png"))
#   fig_tif <- file.path(fig_dir, paste0("Figure2_direct_", trait_set, "_", heatmap_mode, ".tiff"))
  
#   panel_titles <- c("10–16 weeks", "36–38 weeks", "Postpartum")
#   panel_letters <- c("A", "B", "C")
  
#   n_modules <- nrow(moduleTraitCor_ID[[1]])
  
#   fig_width_mm  <- FIG_WIDTH_FULL
#   fig_height_mm <- FIG_HEIGHT_MAX
  
#   draw_figure2 <- function() {
    
#     layout(matrix(1:3, nrow = 3, ncol = 1))
    
#     par(
#       oma = c(0.4, 0.4, 0.4, 0.4),
#       family = "Helvetica"
#     )
    
#     for (set in seq_along(moduleTraitCor_ID)) {
      
#       cor_mat <- moduleTraitCor_ID[[set]]
#       p_mat   <- moduleTraitPvalue_ID[[set]]
#       q_mat   <- moduleTraitQvalue_ID[[set]]
      
#       cor_mat <- cor_mat[, traits_use, drop = FALSE]
#       p_mat   <- p_mat[, traits_use, drop = FALSE]
#       q_mat   <- q_mat[, traits_use, drop = FALSE]
      
#       # Module order
#       if (!is.null(module_order)) {
#         keep_order <- module_order[module_order %in% rownames(cor_mat)]
#         remaining  <- rownames(cor_mat)[!(rownames(cor_mat) %in% keep_order)]
#         final_order <- c(keep_order, remaining)
        
#         cor_mat <- cor_mat[final_order, , drop = FALSE]
#         p_mat   <- p_mat[final_order, , drop = FALSE]
#         q_mat   <- q_mat[final_order, , drop = FALSE]
#       }
      
#       row_labels <- gsub("^ME", "", rownames(cor_mat))
      
#       if (heatmap_mode == "fdr") {
#         sig_mat <- q_mat
#         cutoff <- fdr_cutoff
#         sig_label <- paste0("* BH q < ", cutoff)
#       } else {
#         sig_mat <- p_mat
#         cutoff <- raw_p_cutoff
#         sig_label <- paste0("* p < ", cutoff)
#       }
      
#       textMatrix <- apply(sig_mat, 2, function(x) {
#         sapply(x, function(v) {
#           ifelse(!is.na(v) && v < cutoff, "*", "")
#         })
#       })
      
#       dim(textMatrix) <- dim(cor_mat)
      
#       par(
#         mar = c(4.8, 7.4, 2.2, 4.2),
#         family = "Helvetica"
#       )
      
#       WGCNA::labeledHeatmap(
#         Matrix        = cor_mat,
#         xLabels       = colnames(cor_mat),
#         yLabels       = rownames(cor_mat),
#         ySymbols      = row_labels,
#         colorLabels   = FALSE,
#         colors        = WGCNA::blueWhiteRed(50),
#         textMatrix    = textMatrix,
#         setStdMargins = FALSE,
#         cex.text      = 0.78,
#         textAdj       = c(0.5, 0.8),
#         zlim          = c(-1, 1),
#         main          = paste0(panel_titles[set], " — ", toupper(method)),
#         cex.lab.y     = 0.62,
#         cex.lab.x     = 0.82,
#         plotLegend    = TRUE,
#         legendLabel   = ifelse(method == "pearson", "Pearson r", "Spearman \u03c1")
#       )
      
#       mtext(
#         panel_letters[set],
#         side = 3,
#         line = 0.45,
#         adj = -0.08,
#         cex = 1.2,
#         font = 2
#       )
      
#       mtext(
#         sig_label,
#         side = 1,
#         line = 4.0,
#         adj = 1,
#         cex = 0.65
#       )
#     }
#   }
  

#     convert_pdf_to_rasters <- function(pdf_file, png_file, tiff_file, dpi = 600) {
    
#     message("Converting PDF to high-resolution PNG/TIFF...")
    
#     # Check for pdftocairo
#     pdftocairo_path <- Sys.which("pdftocairo")
    
#     if (pdftocairo_path == "") {
#         warning(
#         "pdftocairo not found. Skipping PNG/TIFF conversion.\n",
#         "Install poppler or use the PDF directly."
#         )
#         return(invisible(NULL))
#     }
    
#     # Temporary output prefix
#     tmp_prefix <- tempfile("figure2_pdf_render_")
    
#     # Render PDF to PNG at high DPI
#     system2(
#         pdftocairo_path,
#         args = c(
#         "-png",
#         "-r", as.character(dpi),
#         pdf_file,
#         tmp_prefix
#         )
#     )
    
#     rendered_png <- paste0(tmp_prefix, "-1.png")
    
#     if (!file.exists(rendered_png)) {
#         warning("PDF-to-PNG conversion failed. Expected file not found: ", rendered_png)
#         return(invisible(NULL))
#     }
    
#     file.copy(rendered_png, png_file, overwrite = TRUE)
    
#     # Convert PNG to TIFF with LZW compression using magick if available
#     if (requireNamespace("magick", quietly = TRUE)) {
#         img <- magick::image_read(rendered_png)
#         magick::image_write(
#         img,
#         path = tiff_file,
#         format = "tiff",
#         compression = "lzw",
#         density = paste0(dpi, "x", dpi)
#         )
#     } else {
#         warning("magick not available. PNG created but TIFF skipped.")
#     }
    
#     message("  PNG:  ", png_file)
#     message("  TIFF: ", tiff_file)
#     }
  
#   # Vector PDF
#   grDevices::cairo_pdf(
#     filename = fig_pdf,
#     width = fig_width_mm / 25.4,
#     height = fig_height_mm / 25.4,
#     family = "Helvetica"
#   )
#   draw_figure2()
#   dev.off()

#     convert_pdf_to_rasters(
#     pdf_file  = fig_pdf,
#     png_file  = fig_png,
#     tiff_file = fig_tif,
#     dpi       = FIGURE_DPI_FINAL
#     )
#     # ── Convert vector PDF to high-resolution PNG/TIFF ───────
#     # This preserves the good PDF layout and gives cleaner raster exports.

#   message("Direct Figure 2 saved:")
#   message("  PDF:  ", fig_pdf)
#   message("  PNG:  ", fig_png)
#   message("  TIFF: ", fig_tif)
  
#   invisible(list(pdf = fig_pdf, png = fig_png, tiff = fig_tif))
# }

# # ── Compute correlations and FDR ──────────────────────────
# message("\nComputing module-trait correlations (", method, ")...")

# trait_family <- make_trait_family_vector(colnames(pheno_data), cfg)

# message("Trait families:")
# print(table(trait_family))

# moduleTraitCor_ID    <- list()
# moduleTraitPvalue_ID <- list()
# moduleTraitQvalue_ID <- list()

# # Load module order after seeing the module names
# first_me <- consensusMEs_ID[[1]]$data
# module_order <- load_module_order(
#   cfg = cfg,
#   wgcna_run_tag = wgcna_run_tag,
#   current_modules = colnames(first_me)
# )

# for (set in 1:All_nSets_ID) {
  
#   me <- consensusMEs_ID[[set]]$data
#   tp <- tp_labels[set]
  
#   ids <- intersect(rownames(me), rownames(pheno_data))
#   message(tp, " | overlap: ", length(ids), " samples")
  
#   if (length(ids) == 0) {
#     stop("No overlapping sample IDs for ", tp)
#   }
  
#   me_sub <- me[ids, , drop = FALSE]
#   pheno_sub <- pheno_data[ids, , drop = FALSE]
  
#   res <- compute_cor(me_sub, pheno_sub, method)
  
#   p_mat <- res$pval
  
#   q_mat <- apply_fdr(
#     p_mat = p_mat,
#     trait_family = trait_family[colnames(p_mat)],
#     method = fdr_method,
#     scope = fdr_scope
#   )
  
#   moduleTraitCor_ID[[set]]    <- res$cor
#   moduleTraitPvalue_ID[[set]] <- p_mat
#   moduleTraitQvalue_ID[[set]] <- q_mat
  
#   message("  FDR q<", fdr_cutoff, ": ", sum(q_mat < fdr_cutoff, na.rm = TRUE),
#           " module-trait cells")
# }

# names(moduleTraitCor_ID)    <- tp_labels
# names(moduleTraitPvalue_ID) <- tp_labels
# names(moduleTraitQvalue_ID) <- tp_labels

# # ── Save all-trait heatmaps and Excel files ───────────────
# heatmap_folders <- character(0)

# for (mode in heatmap_modes) {
  
#   if (mode == "raw") {
#     folder_name <- paste0("heatmaps_raw_p", cutoff_label(raw_p_cutoff))
#   } else if (mode == "fdr") {
#     folder_name <- paste0(
#       "heatmaps_",
#       fdr_scope,
#       "_",
#       fdr_method,
#       "_q",
#       cutoff_label(fdr_cutoff)
#     )
#   } else {
#     stop("Unsupported heatmap_significance mode: ", mode)
#   }
  
#   heatmap_dir <- file.path(out, folder_name)
#   dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)
#   heatmap_folders <- c(heatmap_folders, heatmap_dir)
  
#   for (set in 1:All_nSets_ID) {
    
#     tp  <- tp_labels[set]
#     tpd <- tp_display[set]
    
#     make_heatmap_panel(
#       cor_mat      = moduleTraitCor_ID[[set]],
#       pval_mat     = moduleTraitPvalue_ID[[set]],
#       qval_mat     = moduleTraitQvalue_ID[[set]],
#       title_str    = paste0(tpd, " — ", toupper(method)),
#       pdf_path     = file.path(heatmap_dir, paste0("heatmap_", tp, ".pdf")),
#       fig_panel_id = panel_ids[set],
#       heatmap_mode = mode,
#       module_order = module_order
#     )
#   }
# }

# # Save Excel files with correlations, p-values, q-values
# for (set in 1:All_nSets_ID) {
  
#   tp <- tp_labels[set]
  
#   wb <- openxlsx::createWorkbook()
  
#   openxlsx::addWorksheet(wb, "Correlations")
#   openxlsx::addWorksheet(wb, "P_values")
#   openxlsx::addWorksheet(wb, "BH_Q_values")
  
#   openxlsx::writeData(
#     wb,
#     "Correlations",
#     round(moduleTraitCor_ID[[set]], 3),
#     rowNames = TRUE
#   )
  
#   openxlsx::writeData(
#     wb,
#     "P_values",
#     signif(moduleTraitPvalue_ID[[set]], 3),
#     rowNames = TRUE
#   )
  
#   openxlsx::writeData(
#     wb,
#     "BH_Q_values",
#     signif(moduleTraitQvalue_ID[[set]], 3),
#     rowNames = TRUE
#   )
  
#   openxlsx::saveWorkbook(
#     wb,
#     file.path(out, paste0("ModuleTraitCor_", tp, ".xlsx")),
#     overwrite = TRUE
#   )
# }

# # ── Direct clean Figure 2 for main manuscript ─────────────
# make_direct_figure2_heatmaps(
#   moduleTraitCor_ID    = moduleTraitCor_ID,
#   moduleTraitPvalue_ID = moduleTraitPvalue_ID,
#   moduleTraitQvalue_ID = moduleTraitQvalue_ID,
#   module_order         = module_order,
#   cfg                  = cfg,
#   out_dir              = out
# )

# # ── Save RDS ──────────────────────────────────────────────
# rds_tag <- paste0(wgcna_run_tag, "_", method)

# saveRDS(
#   moduleTraitCor_ID,
#   file.path(cfg$output$processed, paste0("moduleTraitCor_ID_", rds_tag, ".rds"))
# )

# saveRDS(
#   moduleTraitPvalue_ID,
#   file.path(cfg$output$processed, paste0("moduleTraitPvalue_ID_", rds_tag, ".rds"))
# )

# saveRDS(
#   moduleTraitQvalue_ID,
#   file.path(cfg$output$processed, paste0("moduleTraitQvalue_ID_", rds_tag, ".rds"))
# )

# writeLines(
#   paste0(wgcna_run_tag, "|", method),
#   file.path(cfg$output$processed, "current_analysis_tag.txt")
# )

# message("\nScript 08 complete → ", out)
# message("Heatmap folders created:")
# for (hf in unique(heatmap_folders)) {
#   message("  - ", hf)
# }
# message("Figure 2 panels saved to results/figures/panels/fig2/<heatmap_mode>/")
# message("Direct Figure 2 saved to results/figures/final/")
