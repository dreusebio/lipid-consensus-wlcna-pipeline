# R/13_sensitivity_analyses.R
# -------------------------------------------------------
# Sensitivity analysis: Covariate-adjusted module-trait correlations
#
# Reviewer request (R2 comment 4):
#   The parent trial randomized participants to a digital dietary
#   intervention. Any lipid-to-APO or lipid-to-PPWR association may
#   be confounded by intervention assignment and participant characteristics. This script:
#     1. Reports intervention arm/covariate summaries across APO and PPWR groups
#     2. Reruns all module-trait correlations adjusted for config-defined covariates
#        (by residualizing MEs and traits on the configured covariates before correlating)
#     3. Compares covariate-adjusted vs unadjusted findings
#     4. Produces a supplementary covariate-adjusted heatmap figure matching Figure 2 style
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
cor_method <- if (method %in% c("pearson", "spearman", "kendall")) method else "spearman"
if (cor_method != method) {
  message("  Note: adjusted sensitivity uses stats::cor method = ", cor_method,
          " because method = ", method, " is not supported by stats::cor.")
}

# ── Output directory ──────────────────────────────────────
out <- file.path(cfg$output$s14, wgcna_run_tag, method)
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
fdr_method    <- cfg$module_trait$fdr_method   %||% "BH"
fdr_scope     <- cfg$module_trait$fdr_scope    %||% "family"
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

# ── FDR helper functions — same logic as script 08 ─────────
make_trait_family_vector_sens <- function(traits, cfg) {
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

apply_fdr_sens <- function(p_mat, trait_family, method = "BH", scope = "family") {
  q_mat <- matrix(
    NA_real_,
    nrow = nrow(p_mat),
    ncol = ncol(p_mat),
    dimnames = dimnames(p_mat)
  )

  if (scope == "all_traits") {
    q_mat[] <- stats::p.adjust(as.vector(p_mat), method = method)

  } else if (scope == "family") {
    for (fam in unique(trait_family)) {
      fam_traits <- names(trait_family)[trait_family == fam]
      fam_traits <- fam_traits[fam_traits %in% colnames(p_mat)]

      if (length(fam_traits) == 0) next

      p_sub <- p_mat[, fam_traits, drop = FALSE]
      q_sub <- matrix(
        stats::p.adjust(as.vector(p_sub), method = method),
        nrow = nrow(p_sub),
        ncol = ncol(p_sub),
        dimnames = dimnames(p_sub)
      )

      q_mat[, fam_traits] <- q_sub
    }

  } else {
    stop("Unsupported fdr_scope: ", scope,
         ". Use 'family' or 'all_traits'.")
  }

  q_mat
}

trait_family_sens <- make_trait_family_vector_sens(traits_use, cfg)

message("Adjusted sensitivity FDR settings: ",
        "method=", fdr_method,
        " | scope=", fdr_scope,
        " | q cutoff=", fdr_cutoff)


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

# ── Covariate adjustment settings ─────────────────────────
# Config-driven covariates. Add/remove variables in config/config.yml.
#
# Example:
# sensitivity:
#   adjust_covariates:
#     - group
#     - mom_gw_age
#     - bmi_mom
#     - race_ethnicity_binary
#   factor_covariates:
#     - group
#     - race_ethnicity_binary
#   derived_covariates:
#     race_ethnicity_binary:
#       source: race_ethnicity
#       reference_regex: "White"
#       reference_label: "White"
#       other_label: "Non_White"
#
# If this block is absent, the script defaults to intervention arm only
# to preserve the previous behavior.
adjust_covariates <- cfg$sensitivity$adjust_covariates %||%
  cfg$module_trait$adjust_covariates %||%
  c("group")

factor_covariates <- cfg$sensitivity$factor_covariates %||%
  cfg$module_trait$factor_covariates %||%
  character(0)

derived_covariates <- cfg$sensitivity$derived_covariates %||% list()

adjust_covariates <- unique(as.character(unlist(adjust_covariates)))
factor_covariates <- unique(as.character(unlist(factor_covariates)))

if (length(adjust_covariates) == 0) {
  stop("No covariates specified. Add cfg$sensitivity$adjust_covariates in config/config.yml.")
}

# ── Optional derived covariates, useful for sparse race/ethnicity ──
make_derived_covariates <- function(dat, derived_cfg) {
  if (length(derived_cfg) == 0) return(dat)

  for (new_var in names(derived_cfg)) {
    spec <- derived_cfg[[new_var]]
    src  <- spec$source %||% NA_character_
    regex <- spec$reference_regex %||% spec$regex %||% NA_character_
    ref_label <- spec$reference_label %||% "Reference"
    other_label <- spec$other_label %||% "Other"

    if (is.na(src) || !src %in% colnames(dat)) {
      warning("Derived covariate '", new_var, "' skipped: source column not found: ", src)
      next
    }
    if (is.na(regex)) {
      warning("Derived covariate '", new_var, "' skipped: reference_regex not provided.")
      next
    }

    x <- as.character(dat[[src]])
    dat[[new_var]] <- ifelse(
      is.na(x), NA_character_,
      ifelse(grepl(regex, x, ignore.case = TRUE), ref_label, other_label)
    )
  }
  dat
}

demographic_data <- make_derived_covariates(demographic_data, derived_covariates)

missing_covars <- setdiff(adjust_covariates, colnames(demographic_data))
if (length(missing_covars) > 0) {
  stop(
    "The following adjustment covariates are missing from demographic_data: ",
    paste(missing_covars, collapse = ", "),
    "\nAvailable columns include:\n",
    paste(colnames(demographic_data), collapse = ", ")
  )
}

# Coerce configured factor covariates
for (v in intersect(factor_covariates, colnames(demographic_data))) {
  demographic_data[[v]] <- as.factor(demographic_data[[v]])
}

# Also treat character covariates as factors
for (v in adjust_covariates) {
  if (is.character(demographic_data[[v]])) {
    demographic_data[[v]] <- as.factor(demographic_data[[v]])
  }
}

covar_df <- demographic_data[, adjust_covariates, drop = FALSE]
covar_complete_ids <- rownames(covar_df)[stats::complete.cases(covar_df)]

adjustment_label <- paste(adjust_covariates, collapse = " + ")
adjustment_slug  <- paste(gsub("[^A-Za-z0-9]+", "_", adjust_covariates), collapse = "_")
message("Adjustment covariates: ", adjustment_label)
message("Complete covariate data available for n = ", length(covar_complete_ids), " samples")

# Do not include adjustment covariates as traits in the adjusted heatmap.
# Example: if prepregnancy BMI is a covariate, its residual after adjustment
# for itself is not biologically interpretable.
traits_use_before_cov_drop <- traits_use
traits_use <- setdiff(traits_use, adjust_covariates)
dropped_traits <- setdiff(traits_use_before_cov_drop, traits_use)
if (length(dropped_traits) > 0) {
  message("Dropping traits that are also adjustment covariates: ",
          paste(dropped_traits, collapse = ", "))
}

# ── Section 1: Intervention arm distribution table, if arm is available ──
message("\n── Section 1: Covariate summaries across APO and PPWR groups ──")

apo_outcomes  <- c("apo", "apo_hdp", "apo_gdm", "apo_other", "preterm")
ppwr_outcomes <- c("ppwr_e")

wb_cov <- openxlsx::createWorkbook()

if ("group" %in% colnames(demographic_data)) {
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
  if (!is.null(arm_dist_df)) {
    openxlsx::addWorksheet(wb_cov, "Arm_Distribution")
    openxlsx::writeData(wb_cov, "Arm_Distribution", arm_dist_df)
    message("  Arm distribution table added.")
    print(arm_dist_df)
  }
}

covar_summary <- do.call(rbind, lapply(adjust_covariates, function(v) {
  x <- demographic_data[[v]]
  if (is.numeric(x)) {
    data.frame(
      Covariate = v,
      Type = "numeric",
      N_nonmissing = sum(!is.na(x)),
      Summary = sprintf("mean=%.3f; sd=%.3f; median=%.3f; range=%.3f to %.3f",
                        mean(x, na.rm = TRUE), stats::sd(x, na.rm = TRUE),
                        stats::median(x, na.rm = TRUE),
                        min(x, na.rm = TRUE), max(x, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  } else {
    tb <- table(x, useNA = "ifany")
    data.frame(
      Covariate = v,
      Type = "categorical",
      N_nonmissing = sum(!is.na(x)),
      Summary = paste(names(tb), as.integer(tb), sep = "=", collapse = "; "),
      stringsAsFactors = FALSE
    )
  }
}))

openxlsx::addWorksheet(wb_cov, "Covariate_Summary")
openxlsx::writeData(wb_cov, "Covariate_Summary", covar_summary)

openxlsx::saveWorkbook(
  wb_cov,
  file.path(out, "sens_covariate_summary_APO_PPWR.xlsx"),
  overwrite = TRUE
)
message("  Covariate summary workbook saved.")

# ── Section 2: Residualize on configured covariates ───────
residualize_on_covariates <- function(mat, covar_df) {
  ids <- intersect(rownames(mat), rownames(covar_df))
  mat_sub   <- mat[ids, , drop = FALSE]
  covar_sub <- covar_df[ids, , drop = FALSE]

  resid_mat <- apply(mat_sub, 2, function(x) {
    complete <- !is.na(x) & stats::complete.cases(covar_sub)
    out_vec  <- rep(NA_real_, length(x))

    if (sum(complete) > (ncol(covar_sub) + 2)) {
      dat <- data.frame(y = x[complete], covar_sub[complete, , drop = FALSE],
                        check.names = FALSE)
      fit <- tryCatch(stats::lm(y ~ ., data = dat),
                      error = function(e) NULL)
      if (!is.null(fit)) {
        out_vec[complete] <- stats::residuals(fit)
      }
    }
    out_vec
  })

  if (is.null(dim(resid_mat))) {
    resid_mat <- matrix(resid_mat, ncol = 1)
    colnames(resid_mat) <- colnames(mat_sub)
  }

  rownames(resid_mat) <- ids
  colnames(resid_mat) <- colnames(mat_sub)
  resid_mat
}

cor_pvalue_matrix <- function(x, y, method = "spearman") {
  cor_mat <- stats::cor(x, y, method = method, use = "pairwise.complete.obs")
  p_mat   <- matrix(NA_real_, nrow = nrow(cor_mat), ncol = ncol(cor_mat),
                    dimnames = dimnames(cor_mat))
  n_mat   <- matrix(NA_integer_, nrow = nrow(cor_mat), ncol = ncol(cor_mat),
                    dimnames = dimnames(cor_mat))

  for (i in seq_len(ncol(x))) {
    for (j in seq_len(ncol(y))) {
      ok <- stats::complete.cases(x[, i], y[, j])
      n_ij <- sum(ok)
      n_mat[i, j] <- n_ij
      r <- cor_mat[i, j]
      if (!is.na(r) && n_ij > 2) {
        r <- pmin(pmax(r, -0.999999), 0.999999)
        t_stat <- r * sqrt((n_ij - 2) / (1 - r^2))
        p_mat[i, j] <- 2 * stats::pt(-abs(t_stat), df = n_ij - 2)
      }
    }
  }

  list(cor = cor_mat, p = p_mat, n = n_mat)
}

# ── Section 2: Covariate-adjusted correlations ───────────
message("\n── Section 2: Covariate-adjusted module-trait correlations ──")

pheno_numeric <- demographic_data[,
  sapply(demographic_data, is.numeric), drop = FALSE]
pheno_numeric <- pheno_numeric[,
  setdiff(colnames(pheno_numeric), adjust_covariates), drop = FALSE]

adj_cor_ID  <- list()
adj_pval_ID <- list()
adj_qval_ID <- list()
adj_n_ID    <- list()

for (set in seq_len(All_nSets_ID)) {

  tp <- tp_labels[set]
  me <- consensusMEs_ID[[set]]$data

  ids <- Reduce(intersect, list(
    rownames(me),
    rownames(pheno_numeric),
    rownames(covar_df),
    covar_complete_ids
  ))

  message(sprintf("  %s: n=%d samples with complete covariate data", tp, length(ids)))

  me_resid    <- residualize_on_covariates(me[ids, , drop = FALSE],
                                           covar_df[ids, , drop = FALSE])
  trait_resid <- residualize_on_covariates(pheno_numeric[ids, , drop = FALSE],
                                           covar_df[ids, , drop = FALSE])

  traits_in <- intersect(traits_use, colnames(trait_resid))

  cp <- cor_pvalue_matrix(
    me_resid,
    trait_resid[, traits_in, drop = FALSE],
    method = cor_method
  )

  cor_mat  <- cp$cor
  pval_mat <- cp$p
  n_mat    <- cp$n

  qval_mat <- apply_fdr_sens(
    p_mat = pval_mat,
    trait_family = trait_family_sens[colnames(pval_mat)],
    method = fdr_method,
    scope = fdr_scope
  )

  adj_cor_ID[[set]]  <- cor_mat
  adj_pval_ID[[set]] <- pval_mat
  adj_qval_ID[[set]] <- qval_mat
  adj_n_ID[[set]]    <- n_mat

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Correlations")
  openxlsx::addWorksheet(wb, "P_values")
  openxlsx::addWorksheet(wb, "BH_Q_values")
  openxlsx::addWorksheet(wb, "N_pairwise")
  openxlsx::writeData(wb, "Correlations", round(cor_mat,   3), rowNames = TRUE)
  openxlsx::writeData(wb, "P_values",     signif(pval_mat, 3), rowNames = TRUE)
  openxlsx::writeData(wb, "BH_Q_values",  signif(qval_mat, 3), rowNames = TRUE)
  openxlsx::writeData(wb, "N_pairwise",   n_mat, rowNames = TRUE)

  openxlsx::saveWorkbook(
    wb,
    file.path(out, paste0("sens_covariate_adjusted_", tp, ".xlsx")),
    overwrite = TRUE
  )
  message(sprintf("  %s: saved covariate-adjusted correlations", tp))
}

names(adj_cor_ID)  <- tp_labels
names(adj_pval_ID) <- tp_labels
names(adj_qval_ID) <- tp_labels
names(adj_n_ID)    <- tp_labels

saveRDS(adj_cor_ID,  file.path(out, paste0("moduleTraitCor_covariate_adjusted_", adjustment_slug, ".rds")))
saveRDS(adj_pval_ID, file.path(out, paste0("moduleTraitPvalue_covariate_adjusted_", adjustment_slug, ".rds")))
saveRDS(adj_qval_ID, file.path(out, paste0("moduleTraitQvalue_covariate_adjusted_", adjustment_slug, ".rds")))
saveRDS(adj_n_ID,    file.path(out, paste0("moduleTraitN_covariate_adjusted_", adjustment_slug, ".rds")))

# ── Section 3: Comparison unadjusted vs covariate-adjusted ─
message("\n── Section 3: Comparing unadjusted vs covariate-adjusted findings ──")

primary_traits <- c("apo","ppwr_e","apo_hdp","apo_gdm","apo_other","preterm")

comparison <- do.call(rbind, lapply(seq_len(All_nSets_ID), function(set) {
  tp <- tp_labels[set]
  unadj_cor  <- moduleTraitCor_ID[[set]]
  unadj_pval <- moduleTraitPvalue_ID[[set]]
  unadj_qval <- moduleTraitQvalue_ID[[set]]
  adj_cor    <- adj_cor_ID[[set]]
  adj_pval   <- adj_pval_ID[[set]]
  adj_qval   <- adj_qval_ID[[set]]
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
        q_unadj      = signif(unadj_qval[me, tr],  3),
        r_adj        = round(adj_cor[me, tr],      3),
        p_adj        = signif(adj_pval[me, tr],    3),
        q_adj        = signif(adj_qval[me, tr],    3),
        sig_unadj    = unadj_pval[me, tr] < 0.05,
        sig_adj      = adj_pval[me, tr]   < 0.05,
        fdr_unadj    = unadj_qval[me, tr] < fdr_cutoff,
        fdr_adj      = adj_qval[me, tr]   < fdr_cutoff,
        finding_held = (unadj_pval[me, tr] < 0.05) == (adj_pval[me, tr] < 0.05),
        Adjustment   = adjustment_label,
        stringsAsFactors = FALSE
      )
    }))
  }))
}))

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
    All_primary_traits                    = comparison,
    Significant_unadjusted                = comparison_unadj_sig,
    Significant_covariate_adjusted        = comparison_adj_sig,
    Retained_after_covariate_adjust   = comparison_retained,
    Lost_after_covariate_adjustment       = comparison_lost,
    Gained_after_covariate_adjust     = comparison_gained,
    Significant_in_either_analysis        = comparison_sig_either,
    FDR_significant_unadjusted            = comparison[comparison$fdr_unadj & !is.na(comparison$fdr_unadj), ],
    FDR_sig_covariate_adjusted    = comparison[comparison$fdr_adj & !is.na(comparison$fdr_adj), ]
  ),
  file.path(out, "sens_covariate_adjusted_vs_unadjusted_comparison.xlsx"),
  overwrite = TRUE
)

message(sprintf(
  "  Nominally significant associations before covariate adjustment: %d",
  n_unadj_sig))

message(sprintf(
  "  Nominally significant associations after covariate adjustment:  %d",
  n_adj_sig))

message(sprintf(
  "  Original associations retained after covariate adjustment:      %d / %d (%.0f%%)",
  n_retained, n_unadj_sig, pct_retained_of_unadj))

message(sprintf(
  "  Original associations lost after covariate adjustment:          %d / %d (%.0f%%)",
  n_lost, n_unadj_sig, pct_lost_of_unadj))

message(sprintf(
  "  New associations emerging after covariate adjustment:           %d",
  n_gained))

message(sprintf(
  "  Associations significant in either analysis:                    %d",
  n_either))

# ── Section 4: Supplementary heatmap ─────────────────────
message("\n── Section 4: Supplementary heatmap (covariate-adjusted, mirrors Figure 2) ──")

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
# then restricted to modules present in the covariate-adjusted results.

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

# ── Restrict to modules present in covariate-adjusted results ───
current_mods_s13 <- normalize_module_names(rownames(adj_cor_ID[[1]]))
module_order     <- module_order[module_order %in% current_mods_s13]
missing_mods     <- current_mods_s13[!current_mods_s13 %in% module_order]
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

  cor_mat  <- adj_cor_ID[[set_index]]
  pval_mat <- adj_pval_ID[[set_index]]
  qval_mat <- adj_qval_ID[[set_index]]

  traits_in <- intersect(traits_use, colnames(cor_mat))
  cor_mat   <- cor_mat[,  traits_in, drop = FALSE]
  pval_mat  <- pval_mat[, traits_in, drop = FALSE]

  rownames(cor_mat)  <- normalize_module_names(rownames(cor_mat))
  rownames(pval_mat) <- normalize_module_names(rownames(pval_mat))
  rownames(qval_mat) <- normalize_module_names(rownames(qval_mat))

  keep_order <- module_order[module_order %in% rownames(cor_mat)]
  cor_mat    <- cor_mat[keep_order,  , drop = FALSE]
  pval_mat   <- pval_mat[keep_order, , drop = FALSE]
  qval_mat   <- qval_mat[keep_order, , drop = FALSE]

  if (heatmap_mode == "fdr") {
    sig_mat <- qval_mat
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
      name   = ifelse(cor_method=="pearson", "Pearson r", "Spearman \u03c1")) +

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

fig_pdf  <- file.path(fig_dir, "FigureS_covariate_adjusted_heatmap.pdf")
fig_png  <- file.path(fig_dir, "FigureS_covariate_adjusted_heatmap.png")
fig_tiff <- file.path(fig_dir, "FigureS_covariate_adjusted_heatmap.tiff")

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
  "  Before covariate adjustment, %d primary module-trait associations were nominally significant.",
  n_unadj_sig))

message(sprintf(
  "  After covariate adjustment, %d of these original associations remained nominally significant (%.0f%%), while %d were attenuated below p < 0.05.",
  n_retained, pct_retained_of_unadj, n_lost))

if (n_gained > 0) {
  message(sprintf(
    "  In addition, %d associations that were not significant in the unadjusted analysis became nominally significant after covariate adjustment.",
    n_gained))
}

message(sprintf(
  "  Overall, %d associations were nominally significant after covariate adjustment.",
  n_adj_sig))

message(sprintf(
  "  Associations significant in either analysis: %d.",
  n_either))

message("\nScript 14 complete -> ", out)

#This is without FDR 
# # R/13_sensitivity_analyses.R
# # -------------------------------------------------------
# # Sensitivity analysis: Covariate-adjusted module-trait correlations
# #
# # Reviewer request (R2 comment 4):
# #   The parent trial randomized participants to a digital dietary
# #   intervention. Any lipid-to-APO or lipid-to-PPWR association may
# #   be confounded by intervention assignment and participant characteristics. This script:
# #     1. Reports intervention arm/covariate summaries across APO and PPWR groups
# #     2. Reruns all module-trait correlations adjusted for config-defined covariates
# #        (by residualizing MEs and traits on the configured covariates before correlating)
# #     3. Compares covariate-adjusted vs unadjusted findings
# #     4. Produces a supplementary covariate-adjusted heatmap figure matching Figure 2 style
# #
# # Outputs → results/13_sensitivity/<wgcna_run_tag>/<method>/
# # -------------------------------------------------------

# source("R/00_load_packages.R")
# source("R/00_figure_theme.R")

# library(ggplot2)
# library(grid)

# `%||%` <- function(a, b) if (!is.null(a)) a else b

# cfg <- yaml::read_yaml("config/config.yml")

# # ── Read analysis tags ────────────────────────────────────
# wgcna_run_tag    <- trimws(readLines(
#   file.path(cfg$output$processed, "current_wgcna_run_tag.txt")))
# analysis_tag_raw <- trimws(readLines(
#   file.path(cfg$output$processed, "current_analysis_tag.txt")))
# analysis_parts   <- strsplit(analysis_tag_raw, "\\|")[[1]]
# wgcna_run_tag    <- analysis_parts[1]
# method           <- tolower(analysis_parts[2])
# rds_tag          <- paste0(wgcna_run_tag, "_", method)

# message("WGCNA run: ", wgcna_run_tag, " | Method: ", toupper(method))
# cor_method <- if (method %in% c("pearson", "spearman", "kendall")) method else "spearman"
# if (cor_method != method) {
#   message("  Note: adjusted sensitivity uses stats::cor method = ", cor_method,
#           " because method = ", method, " is not supported by stats::cor.")
# }

# # ── Output directory ──────────────────────────────────────
# out <- file.path(cfg$output$s13, wgcna_run_tag, method)
# dir.create(out, showWarnings = FALSE, recursive = TRUE)

# fig_dir <- file.path("results", "figures", "final")
# dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# # ── Load data ─────────────────────────────────────────────
# consensusMEs_ID      <- readRDS(file.path(cfg$output$processed,
#   paste0("consensusMEs_ID_", wgcna_run_tag, ".rds")))
# moduleTraitCor_ID    <- readRDS(file.path(cfg$output$processed,
#   paste0("moduleTraitCor_ID_", rds_tag, ".rds")))
# moduleTraitPvalue_ID <- readRDS(file.path(cfg$output$processed,
#   paste0("moduleTraitPvalue_ID_", rds_tag, ".rds")))
# moduleTraitQvalue_ID <- readRDS(file.path(cfg$output$processed,
#   paste0("moduleTraitQvalue_ID_", rds_tag, ".rds")))
# demographic_data     <- readRDS(file.path(cfg$output$processed,
#   "demographic_data_bmi.rds"))
# exprSize_ID          <- readRDS(file.path(cfg$output$processed,
#   "exprSize_ID.rds"))

# All_nSets_ID <- exprSize_ID$nSets
# tp_labels    <- c("Baseline", "Wk36_38", "Postpartum")
# tp_display   <- c("10\u201316 weeks", "36\u201338 weeks", "Postpartum")

# # Settings — mirror script 08
# raw_p_cutoff  <- cfg$module_trait$raw_p_cutoff %||% 0.05
# fdr_cutoff    <- cfg$module_trait$fdr_cutoff   %||% 0.10
# heatmap_mode  <- cfg$module_trait$figure2_heatmap_mode %||% "raw"
# shared_x_axis <- cfg$module_trait$figure2_shared_x_axis %||% TRUE

# panel_label_case <- tolower(cfg$module_trait$figure2_panel_label_case %||% "lower")
# panel_letters    <- if (panel_label_case == "upper") c("A","B","C") else c("a","b","c")
# panel_titles     <- tp_display

# # Traits — mirror Figure 2 config
# if (cfg$module_trait$figure2_trait_set %||% "primary" == "primary") {
#   traits_use <- unlist(cfg$module_trait$primary_figure2_traits)
# } else {
#   traits_use <- colnames(moduleTraitCor_ID[[1]])
# }
# traits_use <- traits_use[traits_use %in% colnames(moduleTraitCor_ID[[1]])]
# traits_use <- setdiff(traits_use, "group")   # remove arm variable

# message("Traits used: ", paste(traits_use, collapse = ", "))

# # ── Lipid class colors — identical to script 08 ───────────
# LIPID_CLASS_COLORS <- c(
#   "TG"          = "#E63946",
#   "DG"          = "#F4A261",
#   "CE"          = "#8E44AD",
#   "PC/LPC"      = "#1D3D8F",
#   "PE/LPE"      = "#F4A261",
#   "SM"          = "#2A9D8F",
#   "Cer"         = "#E9C46A",
#   "FA"          = "#795548",
#   "Oth"         = "#FFFFFF",
#   "Unknown"     = "#2a756cfb",
#   "Grey module" = "#874788fb"
# )

# LIPID_CLASS_PRIORITY <- c(
#   "TG","DG","CE","PC/LPC","PE/LPE","SM","Cer","FA","Oth","Unknown","Grey module"
# )

# normalize_module_names <- function(x) {
#   x <- as.character(x)
#   x <- gsub("^ME", "", x)
#   paste0("ME", x)
# }

# simplify_lipid_class <- function(x) {
#   x <- gsub("^[0-9]+_", "", as.character(x))
#   x <- gsub("-associated", "", x, ignore.case = TRUE)
#   x <- trimws(x)
#   out <- rep("Unknown", length(x))
#   out[grepl("^TG$|triacyl|triglycer",        x, ignore.case=TRUE)] <- "TG"
#   out[grepl("^DG$|diacyl",                   x, ignore.case=TRUE)] <- "DG"
#   out[grepl("^CE$|cholesteryl",              x, ignore.case=TRUE)] <- "CE"
#   out[grepl("PC/LPC|\\bPC\\b|\\bLPC\\b",    x, ignore.case=TRUE)] <- "PC/LPC"
#   out[grepl("PE/LPE|\\bPE\\b|\\bLPE\\b",    x, ignore.case=TRUE)] <- "PE/LPE"
#   out[grepl("^SM$|sphingomyelin",            x, ignore.case=TRUE)] <- "SM"
#   out[grepl("^Cer$|ceramide",                x, ignore.case=TRUE)] <- "Cer"
#   out[grepl("fatty",                         x, ignore.case=TRUE)] <- "FA"
#   out[grepl("mixed|other",                   x, ignore.case=TRUE)] <- "Oth"
#   out
# }

# # ── Covariate adjustment settings ─────────────────────────
# # Config-driven covariates. Add/remove variables in config/config.yml.
# #
# # Example:
# # sensitivity:
# #   adjust_covariates:
# #     - group
# #     - mom_gw_age
# #     - bmi_mom
# #     - race_ethnicity_binary
# #   factor_covariates:
# #     - group
# #     - race_ethnicity_binary
# #   derived_covariates:
# #     race_ethnicity_binary:
# #       source: race_ethnicity
# #       reference_regex: "White"
# #       reference_label: "White"
# #       other_label: "Non_White"
# #
# # If this block is absent, the script defaults to intervention arm only
# # to preserve the previous behavior.
# adjust_covariates <- cfg$sensitivity$adjust_covariates %||%
#   cfg$module_trait$adjust_covariates %||%
#   c("group")

# factor_covariates <- cfg$sensitivity$factor_covariates %||%
#   cfg$module_trait$factor_covariates %||%
#   character(0)

# derived_covariates <- cfg$sensitivity$derived_covariates %||% list()

# adjust_covariates <- unique(as.character(unlist(adjust_covariates)))
# factor_covariates <- unique(as.character(unlist(factor_covariates)))

# if (length(adjust_covariates) == 0) {
#   stop("No covariates specified. Add cfg$sensitivity$adjust_covariates in config/config.yml.")
# }

# # ── Optional derived covariates, useful for sparse race/ethnicity ──
# make_derived_covariates <- function(dat, derived_cfg) {
#   if (length(derived_cfg) == 0) return(dat)

#   for (new_var in names(derived_cfg)) {
#     spec <- derived_cfg[[new_var]]
#     src  <- spec$source %||% NA_character_
#     regex <- spec$reference_regex %||% spec$regex %||% NA_character_
#     ref_label <- spec$reference_label %||% "Reference"
#     other_label <- spec$other_label %||% "Other"

#     if (is.na(src) || !src %in% colnames(dat)) {
#       warning("Derived covariate '", new_var, "' skipped: source column not found: ", src)
#       next
#     }
#     if (is.na(regex)) {
#       warning("Derived covariate '", new_var, "' skipped: reference_regex not provided.")
#       next
#     }

#     x <- as.character(dat[[src]])
#     dat[[new_var]] <- ifelse(
#       is.na(x), NA_character_,
#       ifelse(grepl(regex, x, ignore.case = TRUE), ref_label, other_label)
#     )
#   }
#   dat
# }

# demographic_data <- make_derived_covariates(demographic_data, derived_covariates)

# missing_covars <- setdiff(adjust_covariates, colnames(demographic_data))
# if (length(missing_covars) > 0) {
#   stop(
#     "The following adjustment covariates are missing from demographic_data: ",
#     paste(missing_covars, collapse = ", "),
#     "\nAvailable columns include:\n",
#     paste(colnames(demographic_data), collapse = ", ")
#   )
# }

# # Coerce configured factor covariates
# for (v in intersect(factor_covariates, colnames(demographic_data))) {
#   demographic_data[[v]] <- as.factor(demographic_data[[v]])
# }

# # Also treat character covariates as factors
# for (v in adjust_covariates) {
#   if (is.character(demographic_data[[v]])) {
#     demographic_data[[v]] <- as.factor(demographic_data[[v]])
#   }
# }

# covar_df <- demographic_data[, adjust_covariates, drop = FALSE]
# covar_complete_ids <- rownames(covar_df)[stats::complete.cases(covar_df)]

# adjustment_label <- paste(adjust_covariates, collapse = " + ")
# adjustment_slug  <- paste(gsub("[^A-Za-z0-9]+", "_", adjust_covariates), collapse = "_")
# message("Adjustment covariates: ", adjustment_label)
# message("Complete covariate data available for n = ", length(covar_complete_ids), " samples")

# # Do not include adjustment covariates as traits in the adjusted heatmap.
# # Example: if prepregnancy BMI is a covariate, its residual after adjustment
# # for itself is not biologically interpretable.
# traits_use_before_cov_drop <- traits_use
# traits_use <- setdiff(traits_use, adjust_covariates)
# dropped_traits <- setdiff(traits_use_before_cov_drop, traits_use)
# if (length(dropped_traits) > 0) {
#   message("Dropping traits that are also adjustment covariates: ",
#           paste(dropped_traits, collapse = ", "))
# }

# # ── Section 1: Intervention arm distribution table, if arm is available ──
# message("\n── Section 1: Covariate summaries across APO and PPWR groups ──")

# apo_outcomes  <- c("apo", "apo_hdp", "apo_gdm", "apo_other", "preterm")
# ppwr_outcomes <- c("ppwr_e")

# wb_cov <- openxlsx::createWorkbook()

# if ("group" %in% colnames(demographic_data)) {
#   arm_dist_list <- lapply(c(apo_outcomes, ppwr_outcomes), function(tr) {
#     if (!tr %in% colnames(demographic_data)) return(NULL)
#     outcome_vec <- demographic_data[[tr]]
#     arm_col     <- demographic_data$group
#     ids <- !is.na(outcome_vec) & !is.na(arm_col)
#     if (sum(ids) == 0) return(NULL)
#     tbl <- table(
#       Arm     = ifelse(arm_col[ids] == 1, "Intervention", "Control"),
#       Outcome = outcome_vec[ids]
#     )
#     df <- as.data.frame(tbl)
#     df$Trait <- tr
#     df
#   })

#   arm_dist_df <- do.call(rbind, Filter(Negate(is.null), arm_dist_list))
#   if (!is.null(arm_dist_df)) {
#     openxlsx::addWorksheet(wb_cov, "Arm_Distribution")
#     openxlsx::writeData(wb_cov, "Arm_Distribution", arm_dist_df)
#     message("  Arm distribution table added.")
#     print(arm_dist_df)
#   }
# }

# covar_summary <- do.call(rbind, lapply(adjust_covariates, function(v) {
#   x <- demographic_data[[v]]
#   if (is.numeric(x)) {
#     data.frame(
#       Covariate = v,
#       Type = "numeric",
#       N_nonmissing = sum(!is.na(x)),
#       Summary = sprintf("mean=%.3f; sd=%.3f; median=%.3f; range=%.3f to %.3f",
#                         mean(x, na.rm = TRUE), stats::sd(x, na.rm = TRUE),
#                         stats::median(x, na.rm = TRUE),
#                         min(x, na.rm = TRUE), max(x, na.rm = TRUE)),
#       stringsAsFactors = FALSE
#     )
#   } else {
#     tb <- table(x, useNA = "ifany")
#     data.frame(
#       Covariate = v,
#       Type = "categorical",
#       N_nonmissing = sum(!is.na(x)),
#       Summary = paste(names(tb), as.integer(tb), sep = "=", collapse = "; "),
#       stringsAsFactors = FALSE
#     )
#   }
# }))

# openxlsx::addWorksheet(wb_cov, "Covariate_Summary")
# openxlsx::writeData(wb_cov, "Covariate_Summary", covar_summary)

# openxlsx::saveWorkbook(
#   wb_cov,
#   file.path(out, "sens_covariate_summary_APO_PPWR.xlsx"),
#   overwrite = TRUE
# )
# message("  Covariate summary workbook saved.")

# # ── Section 2: Residualize on configured covariates ───────
# residualize_on_covariates <- function(mat, covar_df) {
#   ids <- intersect(rownames(mat), rownames(covar_df))
#   mat_sub   <- mat[ids, , drop = FALSE]
#   covar_sub <- covar_df[ids, , drop = FALSE]

#   resid_mat <- apply(mat_sub, 2, function(x) {
#     complete <- !is.na(x) & stats::complete.cases(covar_sub)
#     out_vec  <- rep(NA_real_, length(x))

#     if (sum(complete) > (ncol(covar_sub) + 2)) {
#       dat <- data.frame(y = x[complete], covar_sub[complete, , drop = FALSE],
#                         check.names = FALSE)
#       fit <- tryCatch(stats::lm(y ~ ., data = dat),
#                       error = function(e) NULL)
#       if (!is.null(fit)) {
#         out_vec[complete] <- stats::residuals(fit)
#       }
#     }
#     out_vec
#   })

#   if (is.null(dim(resid_mat))) {
#     resid_mat <- matrix(resid_mat, ncol = 1)
#     colnames(resid_mat) <- colnames(mat_sub)
#   }

#   rownames(resid_mat) <- ids
#   colnames(resid_mat) <- colnames(mat_sub)
#   resid_mat
# }

# cor_pvalue_matrix <- function(x, y, method = "spearman") {
#   cor_mat <- stats::cor(x, y, method = method, use = "pairwise.complete.obs")
#   p_mat   <- matrix(NA_real_, nrow = nrow(cor_mat), ncol = ncol(cor_mat),
#                     dimnames = dimnames(cor_mat))
#   n_mat   <- matrix(NA_integer_, nrow = nrow(cor_mat), ncol = ncol(cor_mat),
#                     dimnames = dimnames(cor_mat))

#   for (i in seq_len(ncol(x))) {
#     for (j in seq_len(ncol(y))) {
#       ok <- stats::complete.cases(x[, i], y[, j])
#       n_ij <- sum(ok)
#       n_mat[i, j] <- n_ij
#       r <- cor_mat[i, j]
#       if (!is.na(r) && n_ij > 2) {
#         r <- pmin(pmax(r, -0.999999), 0.999999)
#         t_stat <- r * sqrt((n_ij - 2) / (1 - r^2))
#         p_mat[i, j] <- 2 * stats::pt(-abs(t_stat), df = n_ij - 2)
#       }
#     }
#   }

#   list(cor = cor_mat, p = p_mat, n = n_mat)
# }

# # ── Section 2: Covariate-adjusted correlations ───────────
# message("\n── Section 2: Covariate-adjusted module-trait correlations ──")

# pheno_numeric <- demographic_data[,
#   sapply(demographic_data, is.numeric), drop = FALSE]
# pheno_numeric <- pheno_numeric[,
#   setdiff(colnames(pheno_numeric), adjust_covariates), drop = FALSE]

# adj_cor_ID  <- list()
# adj_pval_ID <- list()
# adj_n_ID    <- list()

# for (set in seq_len(All_nSets_ID)) {

#   tp <- tp_labels[set]
#   me <- consensusMEs_ID[[set]]$data

#   ids <- Reduce(intersect, list(
#     rownames(me),
#     rownames(pheno_numeric),
#     rownames(covar_df),
#     covar_complete_ids
#   ))

#   message(sprintf("  %s: n=%d samples with complete covariate data", tp, length(ids)))

#   me_resid    <- residualize_on_covariates(me[ids, , drop = FALSE],
#                                            covar_df[ids, , drop = FALSE])
#   trait_resid <- residualize_on_covariates(pheno_numeric[ids, , drop = FALSE],
#                                            covar_df[ids, , drop = FALSE])

#   traits_in <- intersect(traits_use, colnames(trait_resid))

#   cp <- cor_pvalue_matrix(
#     me_resid,
#     trait_resid[, traits_in, drop = FALSE],
#     method = cor_method
#   )

#   cor_mat  <- cp$cor
#   pval_mat <- cp$p
#   n_mat    <- cp$n

#   adj_cor_ID[[set]]  <- cor_mat
#   adj_pval_ID[[set]] <- pval_mat
#   adj_n_ID[[set]]    <- n_mat

#   wb <- openxlsx::createWorkbook()
#   openxlsx::addWorksheet(wb, "Correlations")
#   openxlsx::addWorksheet(wb, "P_values")
#   openxlsx::addWorksheet(wb, "N_pairwise")
#   openxlsx::writeData(wb, "Correlations", round(cor_mat,   3), rowNames = TRUE)
#   openxlsx::writeData(wb, "P_values",     signif(pval_mat, 3), rowNames = TRUE)
#   openxlsx::writeData(wb, "N_pairwise",   n_mat, rowNames = TRUE)

#   openxlsx::saveWorkbook(
#     wb,
#     file.path(out, paste0("sens_covariate_adjusted_", tp, ".xlsx")),
#     overwrite = TRUE
#   )
#   message(sprintf("  %s: saved covariate-adjusted correlations", tp))
# }

# names(adj_cor_ID)  <- tp_labels
# names(adj_pval_ID) <- tp_labels
# names(adj_n_ID)    <- tp_labels

# saveRDS(adj_cor_ID,  file.path(out, paste0("moduleTraitCor_covariate_adjusted_", adjustment_slug, ".rds")))
# saveRDS(adj_pval_ID, file.path(out, paste0("moduleTraitPvalue_covariate_adjusted_", adjustment_slug, ".rds")))
# saveRDS(adj_n_ID,    file.path(out, paste0("moduleTraitN_covariate_adjusted_", adjustment_slug, ".rds")))

# # ── Section 3: Comparison unadjusted vs covariate-adjusted ─
# message("\n── Section 3: Comparing unadjusted vs covariate-adjusted findings ──")

# primary_traits <- c("apo","ppwr_e","apo_hdp","apo_gdm","apo_other","preterm")

# comparison <- do.call(rbind, lapply(seq_len(All_nSets_ID), function(set) {
#   tp <- tp_labels[set]
#   unadj_cor  <- moduleTraitCor_ID[[set]]
#   unadj_pval <- moduleTraitPvalue_ID[[set]]
#   adj_cor    <- adj_cor_ID[[set]]
#   adj_pval   <- adj_pval_ID[[set]]
#   traits_check  <- intersect(intersect(primary_traits, colnames(unadj_cor)),
#                               colnames(adj_cor))
#   modules_check <- intersect(rownames(unadj_cor), rownames(adj_cor))
#   do.call(rbind, lapply(traits_check, function(tr) {
#     do.call(rbind, lapply(modules_check, function(me) {
#       data.frame(
#         Timepoint    = tp,
#         Module       = gsub("^ME", "", me),
#         Trait        = tr,
#         r_unadj      = round(unadj_cor[me, tr],   3),
#         p_unadj      = signif(unadj_pval[me, tr],  3),
#         r_adj        = round(adj_cor[me, tr],      3),
#         p_adj        = signif(adj_pval[me, tr],    3),
#         sig_unadj    = unadj_pval[me, tr] < 0.05,
#         sig_adj      = adj_pval[me, tr]   < 0.05,
#         finding_held = (unadj_pval[me, tr] < 0.05) == (adj_pval[me, tr] < 0.05),
#         Adjustment   = adjustment_label,
#         stringsAsFactors = FALSE
#       )
#     }))
#   }))
# }))

# # ── Define significance categories cleanly ────────────────

# comparison_unadj_sig <- comparison[
#   comparison$sig_unadj & !is.na(comparison$sig_unadj), ]

# comparison_adj_sig <- comparison[
#   comparison$sig_adj & !is.na(comparison$sig_adj), ]

# comparison_retained <- comparison[
#   comparison$sig_unadj & comparison$sig_adj, ]

# comparison_lost <- comparison[
#   comparison$sig_unadj & !comparison$sig_adj, ]

# comparison_gained <- comparison[
#   !comparison$sig_unadj & comparison$sig_adj, ]

# comparison_sig_either <- comparison[
#   comparison$sig_unadj | comparison$sig_adj, ]

# # Counts
# n_unadj_sig <- nrow(comparison_unadj_sig)
# n_adj_sig   <- nrow(comparison_adj_sig)
# n_retained  <- nrow(comparison_retained)
# n_lost      <- nrow(comparison_lost)
# n_gained    <- nrow(comparison_gained)
# n_either    <- nrow(comparison_sig_either)

# pct_retained_of_unadj <- 100 * n_retained / max(n_unadj_sig, 1)
# pct_lost_of_unadj     <- 100 * n_lost     / max(n_unadj_sig, 1)

# # Sanity checks
# stopifnot(n_unadj_sig == n_retained + n_lost)
# stopifnot(n_adj_sig == n_retained + n_gained)

# openxlsx::write.xlsx(
#   list(
#     All_primary_traits                    = comparison,
#     Significant_unadjusted                = comparison_unadj_sig,
#     Significant_covariate_adjusted        = comparison_adj_sig,
#     Retained_after_covariate_adjust   = comparison_retained,
#     Lost_after_covariate_adjustment       = comparison_lost,
#     Gained_after_covariate_adjust     = comparison_gained,
#     Significant_in_either_analysis        = comparison_sig_either
#   ),
#   file.path(out, "sens_covariate_adjusted_vs_unadjusted_comparison.xlsx"),
#   overwrite = TRUE
# )

# message(sprintf(
#   "  Nominally significant associations before covariate adjustment: %d",
#   n_unadj_sig))

# message(sprintf(
#   "  Nominally significant associations after covariate adjustment:  %d",
#   n_adj_sig))

# message(sprintf(
#   "  Original associations retained after covariate adjustment:      %d / %d (%.0f%%)",
#   n_retained, n_unadj_sig, pct_retained_of_unadj))

# message(sprintf(
#   "  Original associations lost after covariate adjustment:          %d / %d (%.0f%%)",
#   n_lost, n_unadj_sig, pct_lost_of_unadj))

# message(sprintf(
#   "  New associations emerging after covariate adjustment:           %d",
#   n_gained))

# message(sprintf(
#   "  Associations significant in either analysis:                    %d",
#   n_either))

# # ── Section 4: Supplementary heatmap ─────────────────────
# message("\n── Section 4: Supplementary heatmap (covariate-adjusted, mirrors Figure 2) ──")

# # ── Helper functions — exact copies from script 08 ───────

# load_module_order_file_s13 <- function(wgcna_run_tag, current_modules) {
#   module_order_file <- file.path(
#     cfg$output$processed,
#     paste0("module_order_by_hub_lipid_class_", wgcna_run_tag, ".rds")
#   )
#   current_modules_norm <- normalize_module_names(current_modules)
#   if (!file.exists(module_order_file)) {
#     message("No hub-lipid-class module order file found. Using current WGCNA module order.")
#     return(current_modules_norm)
#   }
#   message("Loading module order from: ", module_order_file)
#   obj <- readRDS(module_order_file)
#   if (is.data.frame(obj)) {
#     if (!("Module" %in% colnames(obj)))
#       stop("Module order file does not contain a 'Module' column.")
#     if ("PlotOrder" %in% colnames(obj))
#       obj <- obj[order(obj$PlotOrder), , drop = FALSE]
#     module_order <- obj$Module
#   } else {
#     module_order <- obj
#   }
#   module_order <- normalize_module_names(module_order)
#   ordered_keep <- module_order[module_order %in% current_modules_norm]
#   remaining    <- current_modules_norm[!(current_modules_norm %in% ordered_keep)]
#   final_order  <- unique(c(ordered_keep, remaining))
#   message("  Ordered modules loaded: ", length(ordered_keep))
#   message("  Remaining modules appended: ", length(remaining))
#   final_order
# }

# load_module_class_annotation_s13 <- function(wgcna_run_tag, current_modules) {
#   module_class_file <- file.path(
#     cfg$output$processed,
#     paste0("module_order_by_hub_lipid_class_", wgcna_run_tag, ".rds")
#   )
#   current_modules_norm <- normalize_module_names(current_modules)
#   default_df <- data.frame(
#     Module     = current_modules_norm,
#     LipidClass = "Unknown",
#     ClassGroup = "Unknown",
#     ClassColor = unname(LIPID_CLASS_COLORS["Unknown"]),
#     stringsAsFactors = FALSE
#   )
#   if (!file.exists(module_class_file)) {
#     message("No module lipid-class annotation file found. Using white unclassified strip.")
#     return(default_df)
#   }
#   obj <- readRDS(module_class_file)
#   if (!is.data.frame(obj) || !("Module" %in% colnames(obj))) {
#     message("Could not load module class annotation. Using white unclassified strip.")
#     return(default_df)
#   }
#   class_col <- intersect(
#     c("ClassGroup","class_group","DominantClass","dominant_class",
#       "Dominant_Lipid_Class","dominant_lipid_class","LipidClass","lipid_class",
#       "Class","dominant_hub_class","DominantHubClass","dominant_refmet_class"),
#     colnames(obj)
#   )[1]
#   if (is.na(class_col)) {
#     message("Could not find lipid class column. Using white unclassified strip.")
#     return(default_df)
#   }
#   ann <- data.frame(
#     Module     = normalize_module_names(obj$Module),
#     LipidClass = as.character(obj[[class_col]]),
#     stringsAsFactors = FALSE
#   )
#   ann <- ann[!duplicated(ann$Module), , drop = FALSE]
#   ann$ClassGroup <- simplify_lipid_class(ann$LipidClass)
#   ann$ClassColor <- unname(LIPID_CLASS_COLORS[ann$ClassGroup])
#   ann$ClassColor[is.na(ann$ClassColor)] <- unname(LIPID_CLASS_COLORS["Unknown"])
#   out <- merge(default_df[, "Module", drop = FALSE], ann,
#                by = "Module", all.x = TRUE, sort = FALSE)
#   out$LipidClass[is.na(out$LipidClass)] <- "Unknown"
#   out$ClassGroup[is.na(out$ClassGroup)] <- "Unknown"
#   out$ClassColor[is.na(out$ClassColor)] <- unname(LIPID_CLASS_COLORS["Unknown"])
#   # Mark all grey variants as Grey module
#   grey_idx <- grepl("^MEgrey", tolower(out$Module))
#   out$LipidClass[grey_idx] <- "Grey module"
#   out$ClassGroup[grey_idx] <- "Grey module"
#   out$ClassColor[grey_idx] <- unname(LIPID_CLASS_COLORS["Grey module"])
#   message("  Module lipid-class annotations loaded: ", sum(out$ClassGroup != "Unknown"))
#   out
# }

# make_final_module_order_s13 <- function(current_modules,
#                                         module_order_file_order = NULL,
#                                         module_class_df = NULL) {
#   current_modules <- normalize_module_names(current_modules)
#   if (is.null(module_order_file_order)) {
#     module_order_file_order <- current_modules
#   } else {
#     module_order_file_order <- normalize_module_names(module_order_file_order)
#   }
#   module_order_file_order <- module_order_file_order[
#     module_order_file_order %in% current_modules]
#   remaining  <- current_modules[!(current_modules %in% module_order_file_order)]
#   base_order <- unique(c(module_order_file_order, remaining))
#   if (is.null(module_class_df)) return(base_order)
#   class_map  <- setNames(module_class_df$ClassGroup, module_class_df$Module)
#   class_vec  <- class_map[base_order]
#   class_vec[is.na(class_vec)] <- "Unknown"
#   class_rank <- match(class_vec, LIPID_CLASS_PRIORITY)
#   class_rank[is.na(class_rank)] <- length(LIPID_CLASS_PRIORITY)
#   base_order[order(class_rank, seq_along(base_order))]
# }

# # ── Load module order and class annotation using script 08 logic ──
# # Use the full WGCNA module set (not just arm_cor_ID rows) so that
# # the lipid-class priority sorting is computed on the complete set,
# # then restricted to modules present in the covariate-adjusted results.

# all_wgcna_mods <- normalize_module_names(
#   colnames(consensusMEs_ID[[1]]$data)
# )

# module_order_file_order_s13 <- load_module_order_file_s13(
#   wgcna_run_tag   = wgcna_run_tag,
#   current_modules = all_wgcna_mods
# )

# module_class_df <- load_module_class_annotation_s13(
#   wgcna_run_tag   = wgcna_run_tag,
#   current_modules = all_wgcna_mods
# )

# module_order <- make_final_module_order_s13(
#   current_modules         = all_wgcna_mods,
#   module_order_file_order = module_order_file_order_s13,
#   module_class_df         = module_class_df
# )

# # ── Remove WGCNA grey module (exact match only) ───────────
# # grey60 is a legitimate lipid module and must NOT be removed
# grey_idx     <- tolower(gsub("^ME", "", module_order)) == "grey"
# if (any(grey_idx)) {
#   message("  Removing WGCNA grey module from supplementary figure (unassigned lipids)")
# }
# module_order    <- module_order[!grey_idx]
# module_class_df <- module_class_df[
#   tolower(gsub("^ME", "", module_class_df$Module)) != "grey", , drop = FALSE]

# # ── Restrict to modules present in covariate-adjusted results ───
# current_mods_s13 <- normalize_module_names(rownames(adj_cor_ID[[1]]))
# module_order     <- module_order[module_order %in% current_mods_s13]
# missing_mods     <- current_mods_s13[!current_mods_s13 %in% module_order]
# if (length(missing_mods) > 0) {
#   message("  Appending modules not in order file: ",
#           paste(gsub("^ME", "", missing_mods), collapse = ", "))
#   module_order <- unique(c(module_order, missing_mods))
# }

# message("  Final module order: ", length(module_order), " modules")

# # Identical to script 08
# abbrev_map <- c(
#   "TG"="TG","DG"="DG","CE"="CE","PC/LPC"="PC/LPC","PE/LPE"="PE/LPE",
#   "SM"="SM","Cer"="Cer","FA"="FA",
#   "Oth"="Oth","Unknown"="","Grey module"=""
# )
# dark_classes        <- c("TG","PC/LPC","SM","FA")
# # ── Panel plot — exact copy of script 08 make_panel_plot ──
# make_sens_panel <- function(set_index) {

#   show_x_axis <- if (isTRUE(shared_x_axis)) {
#     set_index == All_nSets_ID
#   } else {
#     TRUE
#   }

#   cor_mat  <- adj_cor_ID[[set_index]]
#   pval_mat <- adj_pval_ID[[set_index]]

#   traits_in <- intersect(traits_use, colnames(cor_mat))
#   cor_mat   <- cor_mat[,  traits_in, drop = FALSE]
#   pval_mat  <- pval_mat[, traits_in, drop = FALSE]

#   rownames(cor_mat)  <- normalize_module_names(rownames(cor_mat))
#   rownames(pval_mat) <- normalize_module_names(rownames(pval_mat))

#   keep_order <- module_order[module_order %in% rownames(cor_mat)]
#   cor_mat    <- cor_mat[keep_order,  , drop = FALSE]
#   pval_mat   <- pval_mat[keep_order, , drop = FALSE]

#   if (heatmap_mode == "fdr") {
#     sig_mat <- pval_mat   # adjusted p-values are nominal; use as-is for fdr mode fallback
#     cutoff  <- fdr_cutoff
#   } else {
#     sig_mat <- pval_mat
#     cutoff  <- raw_p_cutoff
#   }

#   module_y <- setNames(rev(seq_along(keep_order)), keep_order)

#   df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
#   colnames(df) <- c("Module", "Trait", "Correlation")
#   sig_df <- as.data.frame(as.table(sig_mat), stringsAsFactors = FALSE)
#   colnames(sig_df) <- c("Module", "Trait", "SigValue")
#   df <- dplyr::left_join(df, sig_df, by = c("Module", "Trait"))

#   df$Module  <- normalize_module_names(df$Module)
#   df$Trait   <- as.character(df$Trait)
#   df$Star    <- ifelse(!is.na(df$SigValue) & df$SigValue < cutoff, "*", "")

#   trait_positions <- seq_along(traits_in)
#   names(trait_positions) <- traits_in
#   df$TraitX  <- trait_positions[df$Trait]
#   df$ModuleY <- module_y[df$Module]

#   # Class annotation
#   ann_df <- module_class_df[, c("Module","ClassGroup","ClassColor"), drop=FALSE]
#   ann_df$Module <- normalize_module_names(ann_df$Module)
#   ann_ordered <- dplyr::left_join(
#     data.frame(Module = keep_order, stringsAsFactors = FALSE),
#     ann_df, by = "Module")
#   ann_ordered$ClassGroup[is.na(ann_ordered$ClassGroup)] <- "Unknown"
#   ann_ordered$ClassColor[is.na(ann_ordered$ClassColor)] <-
#     unname(LIPID_CLASS_COLORS["Unknown"])

#   grey_idx2 <- tolower(gsub("^ME", "", ann_ordered$Module)) == "grey"
#   ann_ordered$ClassGroup[grey_idx2] <- "Grey module"
#   ann_ordered$ClassColor[grey_idx2] <- unname(LIPID_CLASS_COLORS["Grey module"])

#   ann_ordered$ModuleY <- module_y[ann_ordered$Module]

#   r      <- rle(ann_ordered$ClassGroup)
#   ends   <- cumsum(r$lengths)
#   starts <- ends - r$lengths + 1

#   class_blocks <- data.frame(
#     ClassGroup = r$values, start = starts, end = ends,
#     stringsAsFactors = FALSE)
#   class_blocks$ClassColor <- ann_ordered$ClassColor[class_blocks$start]
#   class_blocks$ymin <- pmin(
#     ann_ordered$ModuleY[class_blocks$start],
#     ann_ordered$ModuleY[class_blocks$end]) - 0.5
#   class_blocks$ymax <- pmax(
#     ann_ordered$ModuleY[class_blocks$start],
#     ann_ordered$ModuleY[class_blocks$end]) + 0.5
#   class_blocks$ClassAbbrev <- unname(abbrev_map[class_blocks$ClassGroup])
#   class_blocks$ClassAbbrev[is.na(class_blocks$ClassAbbrev)] <- ""
#   class_blocks$TextColor <- ifelse(
#     class_blocks$ClassGroup %in% dark_classes, "white", "black")

#   p <- ggplot2::ggplot() +

#     ggplot2::geom_rect(
#       data = class_blocks,
#       ggplot2::aes(xmin=0.10, xmax=0.52, ymin=ymin, ymax=ymax),
#       fill=class_blocks$ClassColor, colour="black",
#       linewidth=0.12, inherit.aes=FALSE) +

#     ggplot2::geom_text(
#       data = class_blocks,
#       ggplot2::aes(x=0.31, y=(ymin+ymax)/2,
#                    label=ClassAbbrev, colour=TextColor),
#       angle=90, size=2.3, fontface="bold",
#       family="Helvetica", inherit.aes=FALSE) +

#     ggplot2::scale_colour_identity() +

#     ggplot2::geom_tile(
#       data = df,
#       ggplot2::aes(x=TraitX, y=ModuleY, fill=Correlation),
#       colour="white", linewidth=0.12, width=0.98, height=0.95) +

#     ggplot2::annotate(
#       "text",
#       x      = -6.3,
#       y      = mean(range(df$ModuleY)),
#       label  = panel_titles[set_index],
#       angle  = 90,
#       hjust  = 0.5,
#       vjust  = 1,
#       size   = 5.2,
#       colour = "black",
#       family = "Helvetica") +

#     ggplot2::geom_text(
#       data = df[df$Star != "", , drop=FALSE],
#       ggplot2::aes(x=TraitX, y=ModuleY, label=Star),
#       size=2.8, colour="black",
#       hjust=0.5, vjust=0.5,
#       family="Helvetica") +

#     ggplot2::scale_fill_gradientn(
#       colors = WGCNA::blueWhiteRed(100),
#       limits = c(-1, 1),
#       breaks = c(-1, -0.5, 0, 0.5, 1),
#       name   = ifelse(cor_method=="pearson", "Pearson r", "Spearman \u03c1")) +

#     ggplot2::scale_x_continuous(
#       breaks = trait_positions,
#       labels = if (show_x_axis) traits_in else rep("", length(traits_in)),
#       expand = ggplot2::expansion(mult=c(0.01, 0.02))) +

#     ggplot2::scale_y_continuous(
#       breaks = module_y[keep_order],
#       labels = gsub("^ME", "", keep_order),
#       expand = ggplot2::expansion(mult=c(0.002, 0.002))) +

#     ggplot2::coord_cartesian(
#       xlim = c(0.05, length(traits_in) + 0.5),
#       clip = "off") +

#     ggplot2::labs(
#       title = NULL,
#       tag   = panel_letters[set_index],
#       x = NULL, y = NULL) +

#     ggplot2::theme_minimal(base_family="Helvetica", base_size=8) +

#     ggplot2::theme(
#       plot.title        = ggplot2::element_blank(),
#       plot.tag          = ggplot2::element_text(face="bold", size=13),
#       plot.tag.position = c(-0.085, 1.0),
#       axis.text.x = ggplot2::element_text(
#         angle=45, hjust=1, vjust=1,
#         size  = if (show_x_axis) 8 else 0,
#         colour="black", face="bold"),
#       axis.ticks.x = if (show_x_axis) {
#         ggplot2::element_line(linewidth=0.25)
#       } else {
#         ggplot2::element_blank()
#       },
#       axis.text.y = ggplot2::element_text(
#         size=7, colour="black", face="bold"),
#       panel.grid   = ggplot2::element_blank(),
#       panel.border = ggplot2::element_rect(
#         colour="black", fill=NA, linewidth=0.4),
#       legend.position = "right",
#       legend.title = ggplot2::element_text(size=8, face="bold"),
#       legend.text  = ggplot2::element_text(size=8, face="bold"),
#       plot.margin  = ggplot2::margin(t=4, r=14,
#         b=if (show_x_axis) 8 else 1, l=95)
#     )
#   p
# }

# # ── Draw and save supplementary figure ───────────────────
# p1 <- make_sens_panel(1)
# p2 <- make_sens_panel(2)
# p3 <- make_sens_panel(3)

# draw_supp_figure <- function() {
#   grid::grid.newpage()
#   grid::pushViewport(
#     grid::viewport(
#       layout = grid::grid.layout(
#         nrow    = 3,
#         ncol    = 1,
#         heights = grid::unit(c(1.0, 1.0, 1.1), "null")
#       )
#     )
#   )
#   print(p1, vp = grid::viewport(layout.pos.row=1, layout.pos.col=1))
#   print(p2, vp = grid::viewport(layout.pos.row=2, layout.pos.col=1))
#   print(p3, vp = grid::viewport(layout.pos.row=3, layout.pos.col=1))
#   grid::popViewport()
# }

# fig_pdf  <- file.path(fig_dir, "FigureS_covariate_adjusted_heatmap.pdf")
# fig_png  <- file.path(fig_dir, "FigureS_covariate_adjusted_heatmap.png")
# fig_tiff <- file.path(fig_dir, "FigureS_covariate_adjusted_heatmap.tiff")

# fig_width_mm  <- FIG_WIDTH_FULL
# fig_height_mm <- FIG_HEIGHT_MAX

# grDevices::cairo_pdf(
#   filename = fig_pdf,
#   width    = fig_width_mm  / 25.4,
#   height   = fig_height_mm / 25.4,
#   family   = "Helvetica"
# )
# draw_supp_figure()
# dev.off()

# png(fig_png,
#     width=fig_width_mm, height=fig_height_mm,
#     units="mm", res=FIGURE_DPI_FINAL,
#     type="cairo", family="Helvetica")
# draw_supp_figure()
# dev.off()

# if (requireNamespace("magick", quietly=TRUE)) {
#   img <- magick::image_read(fig_png)
#   magick::image_write(img, path=fig_tiff, format="tiff",
#                       compression="lzw",
#                       density=paste0(FIGURE_DPI_FINAL,"x",FIGURE_DPI_FINAL))
#   message("  TIFF: ", fig_tiff)
# }

# message("Supplementary figure saved:")
# message("  PDF:  ", fig_pdf)
# message("  PNG:  ", fig_png)

# # ── Summary for point-by-point response ──────────────────
# message("\n── Summary for revision response ──")
# message(sprintf(
#   "  Of %d nominally significant module-trait associations (primary outcomes),",
#   n_total))
# message(sprintf(
#   "  %d (%.0f%%) remained nominally significant after adjusting for intervention arm.",
#   n_held, 100 * n_held / max(n_total, 1)))
# if (n_lost > 0)
#   message(sprintf("  %d associations lost significance after adjustment.", n_lost))
# if (n_gained > 0)
#   message(sprintf("  %d additional associations emerged after adjustment.", n_gained))

# message("\nScript 13 complete → ", out)
