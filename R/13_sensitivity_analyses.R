# R/13_sensitivity_analyses.R
# -------------------------------------------------------
# Sensitivity analysis: Arm-adjusted module-trait correlations
#
# Reviewer request (R2 comment 4):
#   The parent trial randomized participants to a digital dietary
#   intervention. Any lipid-to-APO or lipid-to-PPWR association may
#   be confounded by intervention assignment. This script:
#     1. Reports arm distribution across APO and PPWR groups
#     2. Reruns all module-trait correlations adjusted for arm
#        (by residualizing MEs and traits on arm before correlating)
#     3. Compares arm-adjusted vs unadjusted findings
#     4. Produces a supplementary heatmap figure matching Figure 2 style
#
# Outputs → results/13_sensitivity/<wgcna_run_tag>/<method>/
# -------------------------------------------------------

source("R/00_load_packages.R")
source("R/00_figure_theme.R")

library(ggplot2)
library(grid)

`%||%` <- function(a, b) if (!is.null(a)) a else b

cfg <- yaml::read_yaml("config/config.yml")

# ── Read analysis tags ────────────────────────────────────
wgcna_run_tag    <- trimws(readLines(
  file.path(cfg$output$processed, "current_wgcna_run_tag.txt")))
analysis_tag_raw <- trimws(readLines(
  file.path(cfg$output$processed, "current_analysis_tag.txt")))
analysis_parts   <- strsplit(analysis_tag_raw, "\\|")[[1]]
wgcna_run_tag    <- analysis_parts[1]
method           <- tolower(analysis_parts[2])
rds_tag          <- paste0(wgcna_run_tag, "_", method)

message("WGCNA run: ", wgcna_run_tag, " | Method: ", toupper(method))

# ── Output directory ──────────────────────────────────────
out <- file.path(cfg$output$s13, wgcna_run_tag, method)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

fig_dir <- file.path("results", "figures", "final")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load data ─────────────────────────────────────────────
consensusMEs_ID      <- readRDS(file.path(cfg$output$processed,
  paste0("consensusMEs_ID_", wgcna_run_tag, ".rds")))
moduleTraitCor_ID    <- readRDS(file.path(cfg$output$processed,
  paste0("moduleTraitCor_ID_", rds_tag, ".rds")))
moduleTraitPvalue_ID <- readRDS(file.path(cfg$output$processed,
  paste0("moduleTraitPvalue_ID_", rds_tag, ".rds")))
moduleTraitQvalue_ID <- readRDS(file.path(cfg$output$processed,
  paste0("moduleTraitQvalue_ID_", rds_tag, ".rds")))
demographic_data     <- readRDS(file.path(cfg$output$processed,
  "demographic_data_bmi.rds"))
exprSize_ID          <- readRDS(file.path(cfg$output$processed,
  "exprSize_ID.rds"))

All_nSets_ID <- exprSize_ID$nSets
tp_labels    <- c("Baseline", "Wk36_38", "Postpartum")
tp_display   <- c("10\u201316 weeks", "36\u201338 weeks", "Postpartum")

# Settings — mirror script 08
raw_p_cutoff  <- cfg$module_trait$raw_p_cutoff %||% 0.05
fdr_cutoff    <- cfg$module_trait$fdr_cutoff   %||% 0.10
heatmap_mode  <- cfg$module_trait$figure2_heatmap_mode %||% "raw"
shared_x_axis <- cfg$module_trait$figure2_shared_x_axis %||% TRUE

panel_label_case <- tolower(cfg$module_trait$figure2_panel_label_case %||% "lower")
panel_letters    <- if (panel_label_case == "upper") c("A","B","C") else c("a","b","c")
panel_titles     <- tp_display

# Traits — mirror Figure 2 config
if (cfg$module_trait$figure2_trait_set %||% "primary" == "primary") {
  traits_use <- unlist(cfg$module_trait$primary_figure2_traits)
} else {
  traits_use <- colnames(moduleTraitCor_ID[[1]])
}
traits_use <- traits_use[traits_use %in% colnames(moduleTraitCor_ID[[1]])]
traits_use <- setdiff(traits_use, "group")   # remove arm variable

message("Traits used: ", paste(traits_use, collapse = ", "))

# ── Lipid class colors — identical to script 08 ───────────
LIPID_CLASS_COLORS <- c(
  "TG"          = "#E63946",
  "DG"          = "#F4A261",
  "CE"          = "#8E44AD",
  "PC/LPC"      = "#1D3D8F",
  "PE/LPE"      = "#F4A261",
  "SM"          = "#2A9D8F",
  "Cer"         = "#E9C46A",
  "FA"          = "#795548",
  "Oth"         = "#FFFFFF",
  "Unknown"     = "#2a756cfb",
  "Grey module" = "#874788fb"
)

LIPID_CLASS_PRIORITY <- c(
  "TG","DG","CE","PC/LPC","PE/LPE","SM","Cer","FA","Oth","Unknown","Grey module"
)

normalize_module_names <- function(x) {
  x <- as.character(x)
  x <- gsub("^ME", "", x)
  paste0("ME", x)
}

simplify_lipid_class <- function(x) {
  x <- gsub("^[0-9]+_", "", as.character(x))
  x <- gsub("-associated", "", x, ignore.case = TRUE)
  x <- trimws(x)
  out <- rep("Unknown", length(x))
  out[grepl("^TG$|triacyl|triglycer",        x, ignore.case=TRUE)] <- "TG"
  out[grepl("^DG$|diacyl",                   x, ignore.case=TRUE)] <- "DG"
  out[grepl("^CE$|cholesteryl",              x, ignore.case=TRUE)] <- "CE"
  out[grepl("PC/LPC|\\bPC\\b|\\bLPC\\b",    x, ignore.case=TRUE)] <- "PC/LPC"
  out[grepl("PE/LPE|\\bPE\\b|\\bLPE\\b",    x, ignore.case=TRUE)] <- "PE/LPE"
  out[grepl("^SM$|sphingomyelin",            x, ignore.case=TRUE)] <- "SM"
  out[grepl("^Cer$|ceramide",                x, ignore.case=TRUE)] <- "Cer"
  out[grepl("fatty",                         x, ignore.case=TRUE)] <- "FA"
  out[grepl("mixed|other",                   x, ignore.case=TRUE)] <- "Oth"
  out
}

# ── Check arm variable ────────────────────────────────────
if (!"group" %in% colnames(demographic_data)) {
  stop("'group' column not found in demographic_data. Cannot run arm-adjusted sensitivity.")
}

arm_vec <- setNames(demographic_data$group, rownames(demographic_data))

# ── Section 1: Arm distribution table ────────────────────
message("\n── Section 1: Arm distribution across APO and PPWR groups ──")

apo_outcomes  <- c("apo", "apo_hdp", "apo_gdm", "apo_other", "preterm")
ppwr_outcomes <- c("ppwr_e")

arm_dist_list <- lapply(c(apo_outcomes, ppwr_outcomes), function(tr) {
  if (!tr %in% colnames(demographic_data)) return(NULL)
  outcome_vec <- demographic_data[[tr]]
  arm_col     <- demographic_data$group
  ids <- !is.na(outcome_vec) & !is.na(arm_col)
  if (sum(ids) == 0) return(NULL)
  tbl <- table(
    Arm     = ifelse(arm_col[ids] == 1, "Intervention", "Control"),
    Outcome = outcome_vec[ids]
  )
  df <- as.data.frame(tbl)
  df$Trait <- tr
  df
})

arm_dist_df <- do.call(rbind, Filter(Negate(is.null), arm_dist_list))

wb_dist <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_dist, "Arm_Distribution")
openxlsx::writeData(wb_dist, "Arm_Distribution", arm_dist_df)
openxlsx::saveWorkbook(
  wb_dist,
  file.path(out, "sens_arm_distribution_APO_PPWR.xlsx"),
  overwrite = TRUE
)
message("  Arm distribution table saved.")
print(arm_dist_df)

# ── Section 2: Residualize on arm ────────────────────────
residualize_on_arm <- function(mat, arm_vec) {
  ids <- intersect(rownames(mat), names(arm_vec[!is.na(arm_vec)]))
  mat_sub <- mat[ids, , drop = FALSE]
  arm_sub  <- arm_vec[ids]
  resid_mat <- apply(mat_sub, 2, function(x) {
    complete <- !is.na(x)
    out_vec  <- rep(NA_real_, length(x))
    if (sum(complete) > 2) {
      out_vec[complete] <- residuals(lm(x[complete] ~ arm_sub[complete]))
    }
    out_vec
  })
  rownames(resid_mat) <- ids
  resid_mat
}

# ── Section 2: Arm-adjusted correlations ─────────────────
message("\n── Section 2: Arm-adjusted module-trait correlations ──")

pheno_numeric <- demographic_data[,
  sapply(demographic_data, is.numeric), drop = FALSE]
pheno_numeric <- pheno_numeric[,
  setdiff(colnames(pheno_numeric), "group"), drop = FALSE]

arm_cor_ID  <- list()
arm_pval_ID <- list()

for (set in seq_len(All_nSets_ID)) {

  tp <- tp_labels[set]
  me <- consensusMEs_ID[[set]]$data

  ids <- intersect(rownames(me), rownames(pheno_numeric))
  ids <- intersect(ids, names(arm_vec[!is.na(arm_vec)]))

  message(sprintf("  %s: n=%d samples with arm data", tp, length(ids)))

  me_resid    <- residualize_on_arm(me[ids, , drop = FALSE],            arm_vec)
  trait_resid <- residualize_on_arm(pheno_numeric[ids, , drop = FALSE], arm_vec)

  traits_in <- intersect(traits_use, colnames(trait_resid))

  cor_mat <- stats::cor(
    me_resid,
    trait_resid[, traits_in, drop = FALSE],
    method = "spearman",
    use    = "pairwise.complete.obs"
  )

  n        <- nrow(me_resid)
  cor_clip <- pmin(pmax(cor_mat, -0.999999), 0.999999)
  t_stat   <- cor_clip * sqrt((n - 2) / (1 - cor_clip^2))
  pval_mat <- 2 * stats::pt(-abs(t_stat), df = n - 2)
  dimnames(pval_mat) <- dimnames(cor_mat)

  arm_cor_ID[[set]]  <- cor_mat
  arm_pval_ID[[set]] <- pval_mat

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Correlations")
  openxlsx::addWorksheet(wb, "P_values")
  openxlsx::writeData(wb, "Correlations", round(cor_mat,   3), rowNames = TRUE)
  openxlsx::writeData(wb, "P_values",     signif(pval_mat, 3), rowNames = TRUE)
  openxlsx::saveWorkbook(
    wb,
    file.path(out, paste0("sens_arm_adjusted_", tp, ".xlsx")),
    overwrite = TRUE
  )
  message(sprintf("  %s: saved arm-adjusted correlations", tp))
}

names(arm_cor_ID)  <- tp_labels
names(arm_pval_ID) <- tp_labels

# ── Section 3: Comparison unadjusted vs arm-adjusted ─────
message("\n── Section 3: Comparing unadjusted vs arm-adjusted findings ──")

primary_traits <- c("apo","ppwr_e","apo_hdp","apo_gdm","apo_other","preterm")

comparison <- do.call(rbind, lapply(seq_len(All_nSets_ID), function(set) {
  tp <- tp_labels[set]
  unadj_cor  <- moduleTraitCor_ID[[set]]
  unadj_pval <- moduleTraitPvalue_ID[[set]]
  adj_cor    <- arm_cor_ID[[set]]
  adj_pval   <- arm_pval_ID[[set]]
  traits_check  <- intersect(intersect(primary_traits, colnames(unadj_cor)),
                              colnames(adj_cor))
  modules_check <- intersect(rownames(unadj_cor), rownames(adj_cor))
  do.call(rbind, lapply(traits_check, function(tr) {
    do.call(rbind, lapply(modules_check, function(me) {
      data.frame(
        Timepoint    = tp,
        Module       = gsub("^ME", "", me),
        Trait        = tr,
        r_unadj      = round(unadj_cor[me, tr],   3),
        p_unadj      = signif(unadj_pval[me, tr],  3),
        r_adj        = round(adj_cor[me, tr],      3),
        p_adj        = signif(adj_pval[me, tr],    3),
        sig_unadj    = unadj_pval[me, tr] < 0.05,
        sig_adj      = adj_pval[me, tr]   < 0.05,
        finding_held = (unadj_pval[me, tr] < 0.05) == (adj_pval[me, tr] < 0.05),
        stringsAsFactors = FALSE
      )
    }))
  }))
}))

# ── Define significance categories cleanly ────────────────

# ── Define significance categories cleanly ────────────────

comparison_unadj_sig <- comparison[
  comparison$sig_unadj & !is.na(comparison$sig_unadj), ]

comparison_adj_sig <- comparison[
  comparison$sig_adj & !is.na(comparison$sig_adj), ]

comparison_retained <- comparison[
  comparison$sig_unadj & comparison$sig_adj, ]

comparison_lost <- comparison[
  comparison$sig_unadj & !comparison$sig_adj, ]

comparison_gained <- comparison[
  !comparison$sig_unadj & comparison$sig_adj, ]

comparison_sig_either <- comparison[
  comparison$sig_unadj | comparison$sig_adj, ]

# Counts
n_unadj_sig <- nrow(comparison_unadj_sig)
n_adj_sig   <- nrow(comparison_adj_sig)
n_retained  <- nrow(comparison_retained)
n_lost      <- nrow(comparison_lost)
n_gained    <- nrow(comparison_gained)
n_either    <- nrow(comparison_sig_either)

pct_retained_of_unadj <- 100 * n_retained / max(n_unadj_sig, 1)
pct_lost_of_unadj     <- 100 * n_lost     / max(n_unadj_sig, 1)

# Sanity checks
stopifnot(n_unadj_sig == n_retained + n_lost)
stopifnot(n_adj_sig == n_retained + n_gained)

openxlsx::write.xlsx(
  list(
    All_primary_traits              = comparison,
    Significant_unadjusted          = comparison_unadj_sig,
    Significant_arm_adjusted        = comparison_adj_sig,
    Retained_after_arm_adjustment   = comparison_retained,
    Lost_after_arm_adjustment       = comparison_lost,
    Gained_after_arm_adjustment     = comparison_gained,
    Significant_in_either_analysis  = comparison_sig_either
  ),
  file.path(out, "sens_arm_adjusted_vs_unadjusted_comparison.xlsx"),
  overwrite = TRUE
)

message(sprintf(
  "  Nominally significant associations before arm adjustment: %d",
  n_unadj_sig))

message(sprintf(
  "  Nominally significant associations after arm adjustment:  %d",
  n_adj_sig))

message(sprintf(
  "  Original associations retained after arm adjustment:      %d / %d (%.0f%%)",
  n_retained, n_unadj_sig, pct_retained_of_unadj))

message(sprintf(
  "  Original associations lost after arm adjustment:          %d / %d (%.0f%%)",
  n_lost, n_unadj_sig, pct_lost_of_unadj))

message(sprintf(
  "  New associations emerging after arm adjustment:           %d",
  n_gained))

message(sprintf(
  "  Associations significant in either analysis:              %d",
  n_either))

# ── Section 4: Supplementary heatmap ─────────────────────
message("\n── Section 4: Supplementary heatmap (arm-adjusted, mirrors Figure 2) ──")

# ── Helper functions — exact copies from script 08 ───────

load_module_order_file_s13 <- function(wgcna_run_tag, current_modules) {
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
    if (!("Module" %in% colnames(obj)))
      stop("Module order file does not contain a 'Module' column.")
    if ("PlotOrder" %in% colnames(obj))
      obj <- obj[order(obj$PlotOrder), , drop = FALSE]
    module_order <- obj$Module
  } else {
    module_order <- obj
  }
  module_order <- normalize_module_names(module_order)
  ordered_keep <- module_order[module_order %in% current_modules_norm]
  remaining    <- current_modules_norm[!(current_modules_norm %in% ordered_keep)]
  final_order  <- unique(c(ordered_keep, remaining))
  message("  Ordered modules loaded: ", length(ordered_keep))
  message("  Remaining modules appended: ", length(remaining))
  final_order
}

load_module_class_annotation_s13 <- function(wgcna_run_tag, current_modules) {
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
  obj <- readRDS(module_class_file)
  if (!is.data.frame(obj) || !("Module" %in% colnames(obj))) {
    message("Could not load module class annotation. Using white unclassified strip.")
    return(default_df)
  }
  class_col <- intersect(
    c("ClassGroup","class_group","DominantClass","dominant_class",
      "Dominant_Lipid_Class","dominant_lipid_class","LipidClass","lipid_class",
      "Class","dominant_hub_class","DominantHubClass","dominant_refmet_class"),
    colnames(obj)
  )[1]
  if (is.na(class_col)) {
    message("Could not find lipid class column. Using white unclassified strip.")
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
  out <- merge(default_df[, "Module", drop = FALSE], ann,
               by = "Module", all.x = TRUE, sort = FALSE)
  out$LipidClass[is.na(out$LipidClass)] <- "Unknown"
  out$ClassGroup[is.na(out$ClassGroup)] <- "Unknown"
  out$ClassColor[is.na(out$ClassColor)] <- unname(LIPID_CLASS_COLORS["Unknown"])
  # Mark all grey variants as Grey module
  grey_idx <- grepl("^MEgrey", tolower(out$Module))
  out$LipidClass[grey_idx] <- "Grey module"
  out$ClassGroup[grey_idx] <- "Grey module"
  out$ClassColor[grey_idx] <- unname(LIPID_CLASS_COLORS["Grey module"])
  message("  Module lipid-class annotations loaded: ", sum(out$ClassGroup != "Unknown"))
  out
}

make_final_module_order_s13 <- function(current_modules,
                                        module_order_file_order = NULL,
                                        module_class_df = NULL) {
  current_modules <- normalize_module_names(current_modules)
  if (is.null(module_order_file_order)) {
    module_order_file_order <- current_modules
  } else {
    module_order_file_order <- normalize_module_names(module_order_file_order)
  }
  module_order_file_order <- module_order_file_order[
    module_order_file_order %in% current_modules]
  remaining  <- current_modules[!(current_modules %in% module_order_file_order)]
  base_order <- unique(c(module_order_file_order, remaining))
  if (is.null(module_class_df)) return(base_order)
  class_map  <- setNames(module_class_df$ClassGroup, module_class_df$Module)
  class_vec  <- class_map[base_order]
  class_vec[is.na(class_vec)] <- "Unknown"
  class_rank <- match(class_vec, LIPID_CLASS_PRIORITY)
  class_rank[is.na(class_rank)] <- length(LIPID_CLASS_PRIORITY)
  base_order[order(class_rank, seq_along(base_order))]
}

# ── Load module order and class annotation using script 08 logic ──
# Use the full WGCNA module set (not just arm_cor_ID rows) so that
# the lipid-class priority sorting is computed on the complete set,
# then restricted to modules present in the arm-adjusted results.

all_wgcna_mods <- normalize_module_names(
  colnames(consensusMEs_ID[[1]]$data)
)

module_order_file_order_s13 <- load_module_order_file_s13(
  wgcna_run_tag   = wgcna_run_tag,
  current_modules = all_wgcna_mods
)

module_class_df <- load_module_class_annotation_s13(
  wgcna_run_tag   = wgcna_run_tag,
  current_modules = all_wgcna_mods
)

module_order <- make_final_module_order_s13(
  current_modules         = all_wgcna_mods,
  module_order_file_order = module_order_file_order_s13,
  module_class_df         = module_class_df
)

# ── Remove WGCNA grey module (exact match only) ───────────
# grey60 is a legitimate lipid module and must NOT be removed
grey_idx     <- tolower(gsub("^ME", "", module_order)) == "grey"
if (any(grey_idx)) {
  message("  Removing WGCNA grey module from supplementary figure (unassigned lipids)")
}
module_order    <- module_order[!grey_idx]
module_class_df <- module_class_df[
  tolower(gsub("^ME", "", module_class_df$Module)) != "grey", , drop = FALSE]

# ── Restrict to modules present in arm-adjusted results ───
current_mods_s13 <- normalize_module_names(rownames(arm_cor_ID[[1]]))

current_mods_s13 <- current_mods_s13[
  tolower(gsub("^ME", "", current_mods_s13)) != "grey"
]

module_order <- module_order[module_order %in% current_mods_s13]

missing_mods <- current_mods_s13[!current_mods_s13 %in% module_order]
if (length(missing_mods) > 0) {
  message("  Appending modules not in order file: ",
          paste(gsub("^ME", "", missing_mods), collapse = ", "))
  module_order <- unique(c(module_order, missing_mods))
}

message("  Final module order: ", length(module_order), " modules")

# Identical to script 08
abbrev_map <- c(
  "TG"="TG","DG"="DG","CE"="CE","PC/LPC"="PC/LPC","PE/LPE"="PE/LPE",
  "SM"="SM","Cer"="Cer","FA"="FA",
  "Oth"="Oth","Unknown"="","Grey module"=""
)
dark_classes        <- c("TG","PC/LPC","SM","FA")
# ── Panel plot — exact copy of script 08 make_panel_plot ──
make_sens_panel <- function(set_index) {

  show_x_axis <- if (isTRUE(shared_x_axis)) {
    set_index == All_nSets_ID
  } else {
    TRUE
  }

  cor_mat  <- arm_cor_ID[[set_index]]
  pval_mat <- arm_pval_ID[[set_index]]

  traits_in <- intersect(traits_use, colnames(cor_mat))
  cor_mat   <- cor_mat[,  traits_in, drop = FALSE]
  pval_mat  <- pval_mat[, traits_in, drop = FALSE]

  rownames(cor_mat)  <- normalize_module_names(rownames(cor_mat))
  rownames(pval_mat) <- normalize_module_names(rownames(pval_mat))

  keep_order <- module_order[module_order %in% rownames(cor_mat)]
  cor_mat    <- cor_mat[keep_order,  , drop = FALSE]
  pval_mat   <- pval_mat[keep_order, , drop = FALSE]

  if (heatmap_mode == "fdr") {
    sig_mat <- pval_mat   # arm_pval is nominal; use as-is for fdr mode fallback
    cutoff  <- fdr_cutoff
  } else {
    sig_mat <- pval_mat
    cutoff  <- raw_p_cutoff
  }

  module_y <- setNames(rev(seq_along(keep_order)), keep_order)

  df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  colnames(df) <- c("Module", "Trait", "Correlation")
  sig_df <- as.data.frame(as.table(sig_mat), stringsAsFactors = FALSE)
  colnames(sig_df) <- c("Module", "Trait", "SigValue")
  df <- dplyr::left_join(df, sig_df, by = c("Module", "Trait"))

  df$Module  <- normalize_module_names(df$Module)
  df$Trait   <- as.character(df$Trait)
  df$Star    <- ifelse(!is.na(df$SigValue) & df$SigValue < cutoff, "*", "")

  trait_positions <- seq_along(traits_in)
  names(trait_positions) <- traits_in
  df$TraitX  <- trait_positions[df$Trait]
  df$ModuleY <- module_y[df$Module]

  # Class annotation
  ann_df <- module_class_df[, c("Module","ClassGroup","ClassColor"), drop=FALSE]
  ann_df$Module <- normalize_module_names(ann_df$Module)
  ann_ordered <- dplyr::left_join(
    data.frame(Module = keep_order, stringsAsFactors = FALSE),
    ann_df, by = "Module")
  ann_ordered$ClassGroup[is.na(ann_ordered$ClassGroup)] <- "Unknown"
  ann_ordered$ClassColor[is.na(ann_ordered$ClassColor)] <-
    unname(LIPID_CLASS_COLORS["Unknown"])

  grey_idx2 <- tolower(gsub("^ME", "", ann_ordered$Module)) == "grey"
  ann_ordered$ClassGroup[grey_idx2] <- "Grey module"
  ann_ordered$ClassColor[grey_idx2] <- unname(LIPID_CLASS_COLORS["Grey module"])

  ann_ordered$ModuleY <- module_y[ann_ordered$Module]

  r      <- rle(ann_ordered$ClassGroup)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1

  class_blocks <- data.frame(
    ClassGroup = r$values, start = starts, end = ends,
    stringsAsFactors = FALSE)
  class_blocks$ClassColor <- ann_ordered$ClassColor[class_blocks$start]
  class_blocks$ymin <- pmin(
    ann_ordered$ModuleY[class_blocks$start],
    ann_ordered$ModuleY[class_blocks$end]) - 0.5
  class_blocks$ymax <- pmax(
    ann_ordered$ModuleY[class_blocks$start],
    ann_ordered$ModuleY[class_blocks$end]) + 0.5
  class_blocks$ClassAbbrev <- unname(abbrev_map[class_blocks$ClassGroup])
  class_blocks$ClassAbbrev[is.na(class_blocks$ClassAbbrev)] <- ""
  class_blocks$TextColor <- ifelse(
    class_blocks$ClassGroup %in% dark_classes, "white", "black")

  p <- ggplot2::ggplot() +

    ggplot2::geom_rect(
      data = class_blocks,
      ggplot2::aes(xmin=0.10, xmax=0.52, ymin=ymin, ymax=ymax),
      fill=class_blocks$ClassColor, colour="black",
      linewidth=0.12, inherit.aes=FALSE) +

    ggplot2::geom_text(
      data = class_blocks,
      ggplot2::aes(x=0.31, y=(ymin+ymax)/2,
                   label=ClassAbbrev, colour=TextColor),
      angle=90, size=2.3, fontface="bold",
      family="Helvetica", inherit.aes=FALSE) +

    ggplot2::scale_colour_identity() +

    ggplot2::geom_tile(
      data = df,
      ggplot2::aes(x=TraitX, y=ModuleY, fill=Correlation),
      colour="white", linewidth=0.12, width=0.98, height=0.95) +

    ggplot2::annotate(
      "text",
      x      = -6.3,
      y      = mean(range(df$ModuleY)),
      label  = panel_titles[set_index],
      angle  = 90,
      hjust  = 0.5,
      vjust  = 1,
      size   = 5.2,
      colour = "black",
      family = "Helvetica") +

    ggplot2::geom_text(
      data = df[df$Star != "", , drop=FALSE],
      ggplot2::aes(x=TraitX, y=ModuleY, label=Star),
      size=2.8, colour="black",
      hjust=0.5, vjust=0.5,
      family="Helvetica") +

    ggplot2::scale_fill_gradientn(
      colors = WGCNA::blueWhiteRed(100),
      limits = c(-1, 1),
      breaks = c(-1, -0.5, 0, 0.5, 1),
      name   = ifelse(method=="pearson", "Pearson r", "Spearman \u03c1")) +

    ggplot2::scale_x_continuous(
      breaks = trait_positions,
      labels = if (show_x_axis) traits_in else rep("", length(traits_in)),
      expand = ggplot2::expansion(mult=c(0.01, 0.02))) +

    ggplot2::scale_y_continuous(
      breaks = module_y[keep_order],
      labels = gsub("^ME", "", keep_order),
      expand = ggplot2::expansion(mult=c(0.002, 0.002))) +

    ggplot2::coord_cartesian(
      xlim = c(0.05, length(traits_in) + 0.5),
      clip = "off") +

    ggplot2::labs(
      title = NULL,
      tag   = panel_letters[set_index],
      x = NULL, y = NULL) +

    ggplot2::theme_minimal(base_family="Helvetica", base_size=8) +

    ggplot2::theme(
      plot.title        = ggplot2::element_blank(),
      plot.tag          = ggplot2::element_text(face="bold", size=13),
      plot.tag.position = c(-0.085, 1.0),
      axis.text.x = ggplot2::element_text(
        angle=45, hjust=1, vjust=1,
        size  = if (show_x_axis) 8 else 0,
        colour="black", face="bold"),
      axis.ticks.x = if (show_x_axis) {
        ggplot2::element_line(linewidth=0.25)
      } else {
        ggplot2::element_blank()
      },
      axis.text.y = ggplot2::element_text(
        size=7, colour="black", face="bold"),
      panel.grid   = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour="black", fill=NA, linewidth=0.4),
      legend.position = "right",
      legend.title = ggplot2::element_text(size=8, face="bold"),
      legend.text  = ggplot2::element_text(size=8, face="bold"),
      plot.margin  = ggplot2::margin(t=4, r=14,
        b=if (show_x_axis) 8 else 1, l=95)
    )
  p
}

# ── Draw and save supplementary figure ───────────────────
p1 <- make_sens_panel(1)
p2 <- make_sens_panel(2)
p3 <- make_sens_panel(3)

draw_supp_figure <- function() {
  grid::grid.newpage()
  grid::pushViewport(
    grid::viewport(
      layout = grid::grid.layout(
        nrow    = 3,
        ncol    = 1,
        heights = grid::unit(c(1.0, 1.0, 1.1), "null")
      )
    )
  )
  print(p1, vp = grid::viewport(layout.pos.row=1, layout.pos.col=1))
  print(p2, vp = grid::viewport(layout.pos.row=2, layout.pos.col=1))
  print(p3, vp = grid::viewport(layout.pos.row=3, layout.pos.col=1))
  grid::popViewport()
}

fig_pdf  <- file.path(fig_dir, "FigureS_arm_adjusted_heatmap.pdf")
fig_png  <- file.path(fig_dir, "FigureS_arm_adjusted_heatmap.png")
fig_tiff <- file.path(fig_dir, "FigureS_arm_adjusted_heatmap.tiff")

fig_width_mm  <- FIG_WIDTH_FULL
fig_height_mm <- FIG_HEIGHT_MAX

grDevices::cairo_pdf(
  filename = fig_pdf,
  width    = fig_width_mm  / 25.4,
  height   = fig_height_mm / 25.4,
  family   = "Helvetica"
)
draw_supp_figure()
dev.off()

png(fig_png,
    width=fig_width_mm, height=fig_height_mm,
    units="mm", res=FIGURE_DPI_FINAL,
    type="cairo", family="Helvetica")
draw_supp_figure()
dev.off()

if (requireNamespace("magick", quietly=TRUE)) {
  img <- magick::image_read(fig_png)
  magick::image_write(img, path=fig_tiff, format="tiff",
                      compression="lzw",
                      density=paste0(FIGURE_DPI_FINAL,"x",FIGURE_DPI_FINAL))
  message("  TIFF: ", fig_tiff)
}

message("Supplementary figure saved:")
message("  PDF:  ", fig_pdf)
message("  PNG:  ", fig_png)

# ── Summary for point-by-point response ──────────────────
message("\n── Summary for revision response ──")
message(sprintf(
  "  Before adjustment for intervention arm, %d primary module-trait associations were nominally significant.",
  n_unadj_sig))

message(sprintf(
  "  After arm adjustment, %d of these original associations remained nominally significant (%.0f%%), while %d were attenuated below p < 0.05.",
  n_retained, pct_retained_of_unadj, n_lost))

if (n_gained > 0) {
  message(sprintf(
    "  In addition, %d associations that were not significant in the unadjusted analysis became nominally significant after arm adjustment.",
    n_gained))
}

message(sprintf(
  "  Overall, %d associations were nominally significant after arm adjustment.",
  n_adj_sig))
