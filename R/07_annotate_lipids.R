# R/07_annotate_lipids.R
# -------------------------------------------------------
# Annotate hub lipids and all module assignments with:
#   - RefMet classification (Super class, Main class, Sub class)
#   - PubChem CID (via InChiKey from core facility file)
#   - Cross-validation between annotation sources
#
# Three annotation sources:
#   1. Core facility Excel  → identifier, annotation name, InChiKey, m/z, RT, ESI mode
#   2. RefMet results CSV   → lipid class hierarchy (Super/Main/Sub class), HMDB, LM IDs
#   3. PubChem CID CSV      → PubChem CID via InChiKey
#
# Join strategy:
#   annotation name → RefMet (100% match)
#   InChiKey (core) → PubChem CID (96.5% match)
#   NOTE: RefMet PubChem CIDs and direct PubChem CIDs disagree for ~70% of
#   shared records — we flag these and prefer the core InChiKey-based CID
#   as it comes directly from the instrument annotation.
#
# Outputs → results/07_annotations/<run_tag>/
# -------------------------------------------------------

source("R/00_load_packages.R")
source("R/00_figure_theme.R")
cfg <- yaml::read_yaml("config/config.yml")

wgcna_run_tag <- trimws(readLines(
  file.path(cfg$output$processed, "current_wgcna_run_tag.txt")))
message("WGCNA run: ", wgcna_run_tag)

out <- file.path(cfg$output$s07, wgcna_run_tag)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load annotation sources ────────────────────────────
message("Loading annotation sources...")

# Source 1: Core facility Excel — identifier, annotation, InChiKey, m/z, RT, ESI
lipid_excel <- read.xlsx(cfg$data$lipid_excel, sheet = "Data", colNames = FALSE)

# Find header row (contains "identifier")
header_row_idx <- which(apply(lipid_excel, 1, function(r)
  any(grepl("^identifier$", r, ignore.case = TRUE))))[1]

# The header has duplicate "normalized peak height" columns (one per sample)
# and empty strings — keep only the first 7 metadata columns
header     <- as.character(lipid_excel[header_row_idx, ])
data_rows  <- lipid_excel[(header_row_idx + 1):nrow(lipid_excel), 1:7]
colnames(data_rows) <- c("identifier","annotation","ion_species",
                          "InChiKey","mz","ret_time","ESI_mode")
data_rows  <- as.data.frame(data_rows)

# Keep only annotated lipids
core_meta <- data_rows %>%
  filter(!is.na(annotation) & trimws(annotation) != "" &
         trimws(annotation) != "NA") %>%
  mutate(
    annotation = trimws(annotation),
    InChiKey   = trimws(InChiKey),
    mz         = suppressWarnings(as.numeric(mz)),
    ret_time   = suppressWarnings(as.numeric(ret_time))
  ) %>%
  dplyr::select(identifier, annotation, ion_species, InChiKey,
                mz, ret_time, ESI_mode)

message("  Core facility: ", nrow(core_meta), " annotated lipids")

# Source 2: RefMet
refmet <- read.csv(cfg$data$refmet, stringsAsFactors = FALSE) %>%
  mutate(annotation = trimws(annotation))
# Note: read.csv converts spaces to dots in column names automatically
message("  RefMet: ", nrow(refmet), " entries")

# Source 3: PubChem CID via InChiKey
pubchem <- read.csv(cfg$data$pubchem, stringsAsFactors = FALSE) %>%
  mutate(InChiKey = trimws(InChiKey))
message("  PubChem: ", nrow(pubchem), " entries")

# ── 2. Join all sources ───────────────────────────────────
message("\nJoining annotation sources...")

lipid_annot <- core_meta %>%
  # Join RefMet by annotation name (100% match)
  left_join(
    refmet %>%
      dplyr::select(annotation, Standardized.name, Formula, Exact.mass,
                    Super.class, Main.class, Sub.class, Pubchem.CID,
                    CHEBI_ID, HMDB_ID, LM_ID, KEGG_ID, INCHI_KEY, REFMET_ID),
    by = "annotation"
  ) %>%
  # Join PubChem CID via core InChiKey (96.5% match — more reliable)
  left_join(
    pubchem %>% rename(Pubchem_CID_direct = Pubchem_CID),
    by = "InChiKey"
  ) %>%
  # Flag CID disagreements between RefMet and direct PubChem lookup
  mutate(
    Pubchem_CID_refmet = suppressWarnings(as.numeric(Pubchem.CID)),
    CID_agree = case_when(
      is.na(Pubchem_CID_direct) & is.na(Pubchem_CID_refmet) ~ "both_missing",
      is.na(Pubchem_CID_direct)  ~ "only_refmet",
      is.na(Pubchem_CID_refmet)  ~ "only_direct",
      Pubchem_CID_direct == Pubchem_CID_refmet ~ "agree",
      TRUE ~ "disagree"
    ),
    # Use direct (InChiKey-based) CID as primary — more reliable
    Pubchem_CID_final = coalesce(Pubchem_CID_direct, Pubchem_CID_refmet)
  )

message("\nAnnotation coverage:")
message("  RefMet Super class: ",
        sum(!is.na(lipid_annot$Super.class)), "/", nrow(lipid_annot))
message("  RefMet Main class: ",
        sum(!is.na(lipid_annot$Main.class)), "/", nrow(lipid_annot))
message("  PubChem CID (direct): ",
        sum(!is.na(lipid_annot$Pubchem_CID_direct)), "/", nrow(lipid_annot))
message("  PubChem CID (final): ",
        sum(!is.na(lipid_annot$Pubchem_CID_final)), "/", nrow(lipid_annot))
message("\nCID agreement between RefMet and direct lookup:")
print(table(lipid_annot$CID_agree))

# ── 3. Load module assignments and merge ──────────────────
message("\nMerging with module assignments...")

mm_list <- readRDS(file.path(cfg$output$processed,
  paste0("module_membership_", wgcna_run_tag, ".rds")))

# Use baseline module assignments (consistent across timepoints in consensus)
mm_baseline <- mm_list[["Baseline"]] %>%
  dplyr::select(Probe, Module) %>%
  rename(annotation = Probe)

# Join
lipid_full <- lipid_annot %>%
  left_join(mm_baseline, by = "annotation")

message("  Lipids with module assignment: ",
        sum(!is.na(lipid_full$Module)), "/", nrow(lipid_full))

# ── 4. Cross-validation report ────────────────────────────
message("\nGenerating cross-validation report...")

disagree_df <- lipid_full %>%
  filter(CID_agree == "disagree") %>%
  dplyr::select(annotation, InChiKey, Pubchem_CID_direct,
                Pubchem_CID_refmet, Super.class, Main.class)

message("  CID disagreements: ", nrow(disagree_df))
if (nrow(disagree_df) > 0) {
  message("  These lipids have different PubChem CIDs in RefMet vs direct lookup.")
  message("  The direct (InChiKey-based) CID is used as Pubchem_CID_final.")
  message("  Review the validation file to manually verify if needed.")
}

write.csv(disagree_df,
          file.path(out, "cid_disagreements.csv"),
          row.names = FALSE)

# ── 5. Class distribution by module ──────────────────────
message("\nSummarizing lipid class distribution by module...")

class_by_module <- lipid_full %>%
  filter(!is.na(Module) & Module != "grey") %>%
  group_by(Module, Main.class) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(Module, desc(n))

write.csv(class_by_module,
          file.path(out, "lipid_class_by_module.csv"),
          row.names = FALSE)

# ESI mode by module
esi_by_module <- lipid_full %>%
  filter(!is.na(Module) & Module != "grey") %>%
  group_by(Module, ESI_mode) %>%
  summarise(n = n(), .groups = "drop")

write.csv(esi_by_module,
          file.path(out, "esi_mode_by_module.csv"),
          row.names = FALSE)

# ── 6. Save full annotation table ─────────────────────────
final_cols <- c("annotation", "Module", "identifier", "InChiKey",
                "ion_species", "mz", "ret_time", "ESI_mode",
                "Standardized name", "Formula", "Exact mass",
                "Super.class", "Main.class", "Sub.class",
                "Pubchem_CID_final", "Pubchem_CID_direct",
                "Pubchem_CID_refmet", "CID_agree",
                "CHEBI_ID", "HMDB_ID", "LM_ID", "KEGG_ID",
                "INCHI_KEY", "REFMET_ID")

lipid_final <- lipid_full %>%
  dplyr::select(any_of(final_cols)) %>%
  arrange(Module, Super.class, Main.class, annotation)

write.csv(lipid_final,
          file.path(out, "lipid_annotations_full.csv"),
          row.names = FALSE)
write.xlsx(lipid_final,
           file.path(out, "lipid_annotations_full.xlsx"))

# Save as RDS for downstream use
saveRDS(lipid_final,
        file.path(cfg$output$processed,
          paste0("lipid_annotations_", wgcna_run_tag, ".rds")))

# ── 7. Hub lipid annotation tables ────────────────────────
message("\nAnnotating hub lipids...")

tp_labels <- c("Baseline", "Wk36_38", "Postpartum")
for (tp in tp_labels) {
  hub_path <- file.path("results/04_module_membership", wgcna_run_tag,
                        paste0("hub_lipids_top10_", tp, ".csv"))
  if (!file.exists(hub_path)) next

  hubs <- read.csv(hub_path) %>%
    rename(annotation = HubLipid) %>%
    left_join(lipid_final %>%
                dplyr::select(annotation, Super.class, Main.class,
                               Sub.class, Pubchem_CID_final,
                               InChiKey, ESI_mode, mz, HMDB_ID),
              by = "annotation")

  write.xlsx(hubs,
             file.path(out, paste0("hub_lipids_annotated_", tp, ".xlsx")))
  message("  ", tp, ": ", nrow(hubs), " hub lipids annotated")
}


# ── 8. Create hub-based module ordering for heatmaps ──────
# Goal: order modules biologically, with TG-associated modules first.
# Uses annotated top hub lipids when available; falls back to all module lipids.

message("\nCreating hub-based module order...")

# Class group ordering — TGs always first, then DG, CE, PE/LPE, PC/LPC,
# SM, Ceramides, Fatty acyls, Other. Edit the prefix numbers to reorder.
class_group_from_text <- function(x) {
  x <- ifelse(is.na(x), "", x)
  x <- paste(x, collapse = " ")
  dplyr::case_when(
    grepl("Triradylglycerol|triacyl|triglycer|\\bTG\\b", x, ignore.case = TRUE) ~ "01_TG-associated",
    grepl("Diradylglycerol|diacyl|diglycer|\\bDG\\b", x, ignore.case = TRUE) ~ "02_DG-associated",
    grepl("cholest|\\bCE\\b", x, ignore.case = TRUE) ~ "03_Cholesteryl-esters",
    grepl("Glycerophosphoethanolamine|phosphatidylethanolamine|\\bPE\\b|\\bLPE\\b", x, ignore.case = TRUE) ~ "04_PE/LPE",
    grepl("Glycerophosphocholine|phosphatidylcholine|\\bPC\\b|\\bLPC\\b", x, ignore.case = TRUE) ~ "05_PC/LPC",
    grepl("Phosphosphingolipid|sphingomyelin|\\bSM\\b", x, ignore.case = TRUE) ~ "06_Sphingomyelins",
    grepl("Ceramide|\\bCer\\b", x, ignore.case = TRUE) ~ "07_Ceramides",
    grepl("Fatty acid|Fatty ester|Fatty amide|fatty acyl|acylcarnitine", x, ignore.case = TRUE) ~ "08_Fatty-acyls",
    grepl("unknown|other|unclassified", x, ignore.case = TRUE) ~ "99_Other/unknown",
    TRUE ~ "90_Other"
  )
}

hub_annot_all <- list()
for (tp in tp_labels) {
  hub_xlsx <- file.path(out, paste0("hub_lipids_annotated_", tp, ".xlsx"))
  hub_csv  <- file.path(out, paste0("hub_lipids_annotated_", tp, ".csv"))
  hub_raw  <- file.path("results/04_module_membership", wgcna_run_tag,
                        paste0("hub_lipids_top10_", tp, ".csv"))

  if (file.exists(hub_raw)) {
    hubs_tp <- read.csv(hub_raw, stringsAsFactors = FALSE) %>%
      rename(annotation = HubLipid) %>%
      left_join(lipid_final %>%
                  dplyr::select(annotation, Super.class, Main.class, Sub.class),
                by = "annotation") %>%
      mutate(Timepoint = tp)
    hub_annot_all[[tp]] <- hubs_tp
    write.csv(hubs_tp, hub_csv, row.names = FALSE)
  }
}

hub_annot_all <- bind_rows(hub_annot_all)

if (nrow(hub_annot_all) > 0) {
  module_class_summary <- hub_annot_all %>%
    filter(!is.na(Module), Module != "grey") %>%
    mutate(
      ClassText = paste(Super.class, Main.class, Sub.class, sep = " | "),
      ClassGroup = vapply(ClassText, class_group_from_text, character(1))
    ) %>%
    group_by(Module, ClassGroup) %>%
    summarise(
      n_hubs = n(),
      mean_kME = mean(kME, na.rm = TRUE),
      best_rank = min(Rank, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Module, desc(n_hubs), best_rank, desc(mean_kME)) %>%
    group_by(Module) %>%
    slice(1) %>%
    ungroup()
} else {
  message("  No hub annotations found; falling back to all lipids per module.")
  module_class_summary <- lipid_final %>%
    filter(!is.na(Module), Module != "grey") %>%
    mutate(
      ClassText = paste(Super.class, Main.class, Sub.class, sep = " | "),
      ClassGroup = vapply(ClassText, class_group_from_text, character(1))
    ) %>%
    group_by(Module, ClassGroup) %>%
    summarise(n_hubs = n(), mean_kME = NA_real_, best_rank = NA_real_, .groups = "drop") %>%
    arrange(Module, desc(n_hubs)) %>%
    group_by(Module) %>%
    slice(1) %>%
    ungroup()
}

module_order_df <- module_class_summary %>%
  mutate(
    ClassOrder = as.integer(sub("_.*$", "", ClassGroup)),
    Module = as.character(Module)
  ) %>%
  arrange(ClassOrder, ClassGroup, Module) %>%
  mutate(PlotOrder = row_number()) %>%
  dplyr::select(PlotOrder, Module, ClassGroup, n_hubs, mean_kME, best_rank)

write.csv(module_order_df,
          file.path(out, "module_order_by_hub_lipid_class.csv"),
          row.names = FALSE)

saveRDS(module_order_df,
        file.path(cfg$output$processed,
                  paste0("module_order_by_hub_lipid_class_", wgcna_run_tag, ".rds")))

message("  Module order saved: module_order_by_hub_lipid_class.csv")
print(module_order_df)


# ── 9. Lipid class × module heatmaps (Figure 2) ──────────
# heatmap_mode controls what is shown:
#   "count"   → raw number of lipids per class per module
#   "percent" → proportion of each module that is each class (%)
#   "both"    → saves both versions
heatmap_mode <- "both"   # ← change to "count" or "percent" if needed

message("\nGenerating lipid class x module heatmaps (mode: ", heatmap_mode, ")...")

# Build the count matrix: rows = lipid class, cols = module
class_mod_mat <- lipid_full %>%
  filter(!is.na(Module), Module != "grey", !is.na(Main.class)) %>%
  group_by(Main.class, Module) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Module, values_from = n, values_fill = 0)

class_names <- class_mod_mat$Main.class
count_mat   <- as.matrix(class_mod_mat[, -1])
rownames(count_mat) <- class_names

# Order modules by hub lipid class (biological ordering)
if (file.exists(file.path(cfg$output$processed,
    paste0("module_order_by_hub_lipid_class_", wgcna_run_tag, ".rds")))) {
  mod_order <- readRDS(file.path(cfg$output$processed,
    paste0("module_order_by_hub_lipid_class_", wgcna_run_tag, ".rds")))
  ordered_mods <- mod_order$Module[mod_order$Module %in% colnames(count_mat)]
  remaining    <- setdiff(colnames(count_mat), ordered_mods)
  count_mat    <- count_mat[, c(ordered_mods, remaining), drop = FALSE]
}

# Order rows by total count descending
row_totals <- rowSums(count_mat)
count_mat  <- count_mat[order(row_totals, decreasing = TRUE), ]

# Percent matrix
pct_mat <- sweep(count_mat, 2, colSums(count_mat), "/") * 100

# WGCNA module colors — used as column color bars at top of heatmap
# R's built-in color names match WGCNA module names directly
module_color_bar <- function(module_names) {
  # Most WGCNA color names are valid R colors; handle exceptions
  color_map <- setNames(module_names, module_names)
  color_map["grey60"]       <- "grey60"
  color_map["darkturquoise"]<- "darkturquoise"
  color_map["darkolivegreen"]<- "darkolivegreen"
  color_map["saddlebrown"]  <- "saddlebrown"
  color_map["paleturquoise"]<- "paleturquoise"
  color_map["yellowgreen"]  <- "yellowgreen"
  color_map["steelblue"]    <- "steelblue"
  color_map["darkmagenta"]  <- "darkmagenta"
  color_map["darkorange"]   <- "darkorange"
  color_map["darkred"]      <- "darkred"
  color_map["darkgreen"]    <- "darkgreen"
  color_map["darkgrey"]     <- "darkgrey"
  color_map["sienna3"]      <- "sienna3"
  color_map["violet"]       <- "violet"
  color_map["white"]        <- "grey90"   # white module → light grey so visible
  color_map
}

# ── Display name mapping (figures only — data tables unchanged) ───
# Maps RefMet Main.class names to reader-friendly labels for plots
# Excel/CSV outputs always use original RefMet names
LIPID_CLASS_DISPLAY <- c(
  # Glycerolipids
  "Triradylglycerols"                        = "Triglycerides (TG)",
  "Diradylglycerols"                         = "Diglycerides (DG)",
  "Monoradylglycerols"                       = "Monoglycerides (MG)",
  "Glycosyldiradylglycerols"                 = "Glycosyldiacylglycerols (GDG)",
  # Glycerophospholipids
  "Glycerophosphocholines"                   = "Phosphatidylcholines (PC/LPC)",
  "Glycerophosphoethanolamines"              = "Phosphatidylethanolamines (PE/LPE)",
  "Glycerophosphoinositols"                  = "Phosphatidylinositols (PI/LPI)",
  "Glycerophosphoserines"                    = "Phosphatidylserines (PS)",
  "Glycerophosphatidylethanol"               = "Phosphatidylethanol (PEtOH)",
  "Glycerophosphates"                        = "Glycerophosphates (GP)",
  # Sphingolipids
  "Phosphosphingolipids"                     = "Sphingomyelins (SM)",
  "Neutral_glycosphingolipids"               = "Neutral Glycosphingolipids (nGSL)",
  "Glycosphingolipids"                       = "Glycosphingolipids (GSL)",
  "Ceramides"                                = "Ceramides (Cer)",
  "Sphingoid bases"                          = "Sphingoid Bases (SPB)",
  # Fatty Acyls
  "Fatty acids"                              = "Fatty Acids (FA)",
  "Fatty esters"                             = "Fatty Esters (FE)",
  "Fatty amides"                             = "Fatty Amides (FAM)",
  "Eicosanoids"                              = "Eicosanoids (EIC)",
  "Docosanoids"                              = "Docosanoids (DOC)",
  "Octadecanoids"                            = "Octadecanoids (OCT)",
  "Keto acids"                               = "Keto Acids (KA)",
  # Sterol Lipids
  "Sterol esters"                            = "Sterol Esters (SE)",
  "Sterols"                                  = "Sterols (ST)",
  # Other
  "Perfluoroalkyl_Phosphate_Acid_derivative" = "PFAS Derivatives"
)

# Helper: apply display names to a vector of class names
apply_display_names <- function(x) {
  ifelse(x %in% names(LIPID_CLASS_DISPLAY),
         LIPID_CLASS_DISPLAY[x],
         x)   # fallback: keep original if not in map
}

make_class_heatmap <- function(mat, title_str, legend_label,
                                fmt_fn = function(x) ifelse(x == 0, "", as.character(x))) {
  max_val    <- max(mat, na.rm = TRUE)
  n_rows     <- nrow(mat)
  n_cols     <- ncol(mat)
  mod_names  <- colnames(mat)
  mod_colors <- module_color_bar(mod_names)

  # ── Main tile data ────────────────────────────────────────
  # Apply display names to row labels (figures only)
  display_rownames <- apply_display_names(rownames(mat))
  mat_display <- mat
  rownames(mat_display) <- display_rownames

  df <- as.data.frame(mat_display) %>%
    tibble::rownames_to_column("LipidClass") %>%
    tidyr::pivot_longer(-LipidClass, names_to = "Module", values_to = "Value") %>%
    mutate(
      LipidClass  = factor(LipidClass, levels = rev(display_rownames)),
      Module      = factor(Module, levels = mod_names),
      Label       = fmt_fn(Value),
      TextColor   = ifelse(Value > max_val * 0.55, "white", "black"),
      # Zebra stripe: alternating row backgrounds
      RowIdx      = as.integer(LipidClass),
      RowStripe   = ifelse(RowIdx %% 2 == 0, "#F7F7F7", "#FFFFFF")
    )

  # ── Module color bar data (row above x-axis labels) ───────
  color_bar_df <- data.frame(
    Module    = factor(mod_names, levels = mod_names),
    FillColor = unname(mod_colors[mod_names]),
    y         = "Module color"
  )

  # ── Main plot ─────────────────────────────────────────────
  p <- ggplot(df, aes(x = Module, y = LipidClass)) +
    # Zebra stripe background (full row width)
    geom_tile(aes(fill = RowStripe), alpha = 0.6, colour = NA) +
    scale_fill_identity() +
    ggnewscale::new_scale_fill() +
    # Value tiles
    geom_tile(aes(fill = Value), colour = NA,
              width = 0.95, height = 0.95) +
    # Vertical grid lines between every module column
    geom_vline(
      xintercept = seq(0.5, n_cols + 0.5, by = 1),
      colour = "grey80", linewidth = 0.25
    ) +
    # Horizontal grid lines between every lipid class row
    geom_hline(
      yintercept = seq(0.5, n_rows + 0.5, by = 1),
      colour = "grey80", linewidth = 0.25
    ) +
    # Cell text
    geom_text(aes(label = Label, colour = TextColor),
              size = FS_SMALL / ggplot2::.pt, fontface = "bold") +
    scale_colour_identity() +
    scale_fill_gradientn(
      colours  = c("#FFFFFF", "#FDE8D8", "#F4A261", "#E63946", "#8B0000"),
      values   = scales::rescale(c(0, 0.2, 0.5, 0.75, 1)),
      na.value = "grey97",
      name     = legend_label
    ) +
    # Outer border
    annotate("rect",
             xmin = 0.5, xmax = n_cols + 0.5,
             ymin = 0.5, ymax = n_rows + 0.5,
             fill = NA, colour = "grey50", linewidth = 0.4) +
    labs(title = title_str, x = NULL, y = "Lipid class") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_growell() +
    theme(
      axis.text.x     = element_text(angle = 50, hjust = 1,
                                     size  = FS_SMALL - 0.5, colour = "black"),
      axis.text.y     = element_text(size = FS_SMALL, colour = "black"),
      axis.title.y    = element_text(size = FS_BASE),
      axis.ticks      = element_blank(),
      legend.position = "right",
      panel.border    = element_blank(),
      panel.grid      = element_blank()
    )

  # ── Module color bar aligned to heatmap columns ───────────
  p_bar <- ggplot(color_bar_df, aes(x = Module, y = y, fill = FillColor)) +
    geom_tile(colour = "white", linewidth = 0.5, height = 0.85) +
    scale_fill_identity() +
    scale_x_discrete(limits = mod_names, expand = c(0, 0)) +
    labs(x = "Module", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      axis.text.x  = element_blank(),
      axis.title.x = element_text(size = FS_BASE, vjust = -0.5),
      plot.margin  = unit(c(2, 0, 0, 0), "pt")
    )

  # Stack: color bar on top of heatmap, aligned
  p_bar / p +
    patchwork::plot_layout(heights = c(0.05, 1)) &
    theme(plot.margin = unit(c(1, 1, 1, 1), "pt"))
}

# Save based on mode
save_heatmap <- function(p, suffix, panel_id) {
  fname <- file.path(out, paste0("heatmap_lipid_class_module_", suffix, ".pdf"))
  ggsave(fname, plot = p,
         width  = FIG_WIDTH_FULL, height = 120, units = "mm",
         device = cairo_pdf)
  ggsave(gsub(".pdf", ".png", fname), plot = p,
         width = FIG_WIDTH_FULL, height = 120, units = "mm",
         dpi   = FIGURE_DPI)
  save_panel(p, "fig2", panel_id,
             width_mm = FIG_WIDTH_FULL, height_mm = 120)
  message("  Saved: heatmap_lipid_class_module_", suffix, ".pdf → fig2_", panel_id)
}

if (heatmap_mode %in% c("count", "both")) {
  p_count <- make_class_heatmap(
    count_mat,
    title_str    = "Lipid class distribution by module (count)",
    legend_label = "No. lipids",
    fmt_fn       = function(x) ifelse(x == 0, "", as.character(x))
  )
  save_heatmap(p_count, "count", if (heatmap_mode == "both") "D" else "D")
}

if (heatmap_mode %in% c("percent", "both")) {
  p_pct <- make_class_heatmap(
    pct_mat,
    title_str    = "Lipid class distribution by module (%)",
    legend_label = "% of module",
    fmt_fn       = function(x) ifelse(x < 1, "", paste0(round(x), "%"))
  )
  save_heatmap(p_pct, "percent", if (heatmap_mode == "both") "E" else "D")
}

# ── 10. Annotated lipid overview bar chart (Figure 2A) ────
# Bar chart showing number of identified lipids per Main class
# coloured by Super class — overview of the full lipidome detected
message("\nGenerating annotated lipid overview chart...")

# Colorblind-safe palette for lipid super classes
superclass_colors <- c(
  "Fatty Acyls"            = "#E63946",
  "Glycerophospholipids"   = "#1D3D8F",
  "Sphingolipids"          = "#F4A261",
  "Sterol Lipids"          = "#2A9D8F",
  "Glycerolipids"          = "#E9C46A",
  "Prenol Lipids"          = "#264653",
  "Saccharolipids"         = "#A8DADC",
  "Polyketides"            = "#457B9D",
  "Other"                  = "#CCCCCC"
)

lipid_overview <- lipid_full %>%
  mutate(
    Super.class = ifelse(is.na(Super.class) | Super.class == "",
                          "Other", Super.class),
    Main.class  = ifelse(is.na(Main.class) | Main.class == "",
                          "Unknown", Main.class)
  ) %>%
  group_by(Super.class, Main.class) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(Super.class, desc(n)) %>%
  mutate(
    # Apply display names before making factor so levels are correct
    Main.class_display = apply_display_names(as.character(Main.class)),
    Main.class_display = factor(Main.class_display,
                                 levels = rev(unique(Main.class_display)))
  )

p_overview <- ggplot(lipid_overview,
                      aes(x = n, y = Main.class_display, fill = Super.class)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = superclass_colors,
                    name   = "Lipid super class") +
  labs(title = "Identified lipids by class",
       x = "Number of lipids", y = NULL) +
  theme_growell(grid = "y") +
  theme(legend.position = "right",
        legend.text     = element_text(size = FS_SMALL),
        axis.text.y     = element_text(size = FS_SMALL))

# Save PDF (using pdf device with Helvetica — cairo_pdf may not be available)
pdf(file.path(out, "annotated_lipid_overview.pdf"),
    width  = FIG_WIDTH_FULL / 25.4,
    height = 110 / 25.4,
    family = "Helvetica")
print(p_overview)
dev.off()

# Also save PNG for easy viewing
ggsave(file.path(out, "annotated_lipid_overview.png"),
       plot  = p_overview,
       width = FIG_WIDTH_FULL, height = 110, units = "mm",
       dpi   = FIGURE_DPI)

save_panel(p_overview, "fig2", "A",
           width_mm = FIG_WIDTH_FULL, height_mm = 100)
message("  Saved: annotated_lipid_overview.pdf → fig2_A")

message("\nScript 07 complete → ", out)
message("  Full annotation: lipid_annotations_full.xlsx")
message("  CID disagreements to review: cid_disagreements.csv (", nrow(disagree_df), " lipids)")
# # R/07_annotate_lipids.R
# # -------------------------------------------------------
# # Annotate hub lipids and all module assignments with:
# #   - RefMet classification (Super class, Main class, Sub class)
# #   - PubChem CID (via InChiKey from core facility file)
# #   - Cross-validation between annotation sources
# #
# # Three annotation sources:
# #   1. Core facility Excel  → identifier, annotation name, InChiKey, m/z, RT, ESI mode
# #   2. RefMet results CSV   → lipid class hierarchy (Super/Main/Sub class), HMDB, LM IDs
# #   3. PubChem CID CSV      → PubChem CID via InChiKey
# #
# # Join strategy:
# #   annotation name → RefMet (100% match)
# #   InChiKey (core) → PubChem CID (96.5% match)
# #   NOTE: RefMet PubChem CIDs and direct PubChem CIDs disagree for ~70% of
# #   shared records — we flag these and prefer the core InChiKey-based CID
# #   as it comes directly from the instrument annotation.
# #
# # Outputs → results/07_annotations/<run_tag>/
# # -------------------------------------------------------

# source("R/00_load_packages.R")
# source("R/00_figure_theme.R")
# cfg <- yaml::read_yaml("config/config.yml")

# wgcna_run_tag <- trimws(readLines(
#   file.path(cfg$output$processed, "current_wgcna_run_tag.txt")))
# message("WGCNA run: ", wgcna_run_tag)

# out <- file.path(cfg$output$s07, wgcna_run_tag)
# dir.create(out, showWarnings = FALSE, recursive = TRUE)

# # ── 1. Load annotation sources ────────────────────────────
# message("Loading annotation sources...")

# # Source 1: Core facility Excel — identifier, annotation, InChiKey, m/z, RT, ESI
# lipid_excel <- read.xlsx(cfg$data$lipid_excel, sheet = "Data", colNames = FALSE)

# # Find header row (contains "identifier")
# header_row_idx <- which(apply(lipid_excel, 1, function(r)
#   any(grepl("^identifier$", r, ignore.case = TRUE))))[1]

# # The header has duplicate "normalized peak height" columns (one per sample)
# # and empty strings — keep only the first 7 metadata columns
# header     <- as.character(lipid_excel[header_row_idx, ])
# data_rows  <- lipid_excel[(header_row_idx + 1):nrow(lipid_excel), 1:7]
# colnames(data_rows) <- c("identifier","annotation","ion_species",
#                           "InChiKey","mz","ret_time","ESI_mode")
# data_rows  <- as.data.frame(data_rows)

# # Keep only annotated lipids
# core_meta <- data_rows %>%
#   filter(!is.na(annotation) & trimws(annotation) != "" &
#          trimws(annotation) != "NA") %>%
#   mutate(
#     annotation = trimws(annotation),
#     InChiKey   = trimws(InChiKey),
#     mz         = suppressWarnings(as.numeric(mz)),
#     ret_time   = suppressWarnings(as.numeric(ret_time))
#   ) %>%
#   dplyr::select(identifier, annotation, ion_species, InChiKey,
#                 mz, ret_time, ESI_mode)

# message("  Core facility: ", nrow(core_meta), " annotated lipids")

# # Source 2: RefMet
# refmet <- read.csv(cfg$data$refmet, stringsAsFactors = FALSE) %>%
#   mutate(annotation = trimws(annotation))
# # Note: read.csv converts spaces to dots in column names automatically
# message("  RefMet: ", nrow(refmet), " entries")

# # Source 3: PubChem CID via InChiKey
# pubchem <- read.csv(cfg$data$pubchem, stringsAsFactors = FALSE) %>%
#   mutate(InChiKey = trimws(InChiKey))
# message("  PubChem: ", nrow(pubchem), " entries")

# # ── 2. Join all sources ───────────────────────────────────
# message("\nJoining annotation sources...")

# lipid_annot <- core_meta %>%
#   # Join RefMet by annotation name (100% match)
#   left_join(
#     refmet %>%
#       dplyr::select(annotation, Standardized.name, Formula, Exact.mass,
#                     Super.class, Main.class, Sub.class, Pubchem.CID,
#                     CHEBI_ID, HMDB_ID, LM_ID, KEGG_ID, INCHI_KEY, REFMET_ID),
#     by = "annotation"
#   ) %>%
#   # Join PubChem CID via core InChiKey (96.5% match — more reliable)
#   left_join(
#     pubchem %>% rename(Pubchem_CID_direct = Pubchem_CID),
#     by = "InChiKey"
#   ) %>%
#   # Flag CID disagreements between RefMet and direct PubChem lookup
#   mutate(
#     Pubchem_CID_refmet = suppressWarnings(as.numeric(Pubchem.CID)),
#     CID_agree = case_when(
#       is.na(Pubchem_CID_direct) & is.na(Pubchem_CID_refmet) ~ "both_missing",
#       is.na(Pubchem_CID_direct)  ~ "only_refmet",
#       is.na(Pubchem_CID_refmet)  ~ "only_direct",
#       Pubchem_CID_direct == Pubchem_CID_refmet ~ "agree",
#       TRUE ~ "disagree"
#     ),
#     # Use direct (InChiKey-based) CID as primary — more reliable
#     Pubchem_CID_final = coalesce(Pubchem_CID_direct, Pubchem_CID_refmet)
#   )

# message("\nAnnotation coverage:")
# message("  RefMet Super class: ",
#         sum(!is.na(lipid_annot$Super.class)), "/", nrow(lipid_annot))
# message("  RefMet Main class: ",
#         sum(!is.na(lipid_annot$Main.class)), "/", nrow(lipid_annot))
# message("  PubChem CID (direct): ",
#         sum(!is.na(lipid_annot$Pubchem_CID_direct)), "/", nrow(lipid_annot))
# message("  PubChem CID (final): ",
#         sum(!is.na(lipid_annot$Pubchem_CID_final)), "/", nrow(lipid_annot))
# message("\nCID agreement between RefMet and direct lookup:")
# print(table(lipid_annot$CID_agree))

# # ── 3. Load module assignments and merge ──────────────────
# message("\nMerging with module assignments...")

# mm_list <- readRDS(file.path(cfg$output$processed,
#   paste0("module_membership_", wgcna_run_tag, ".rds")))

# # Use baseline module assignments (consistent across timepoints in consensus)
# mm_baseline <- mm_list[["Baseline"]] %>%
#   dplyr::select(Probe, Module) %>%
#   rename(annotation = Probe)

# # Join
# lipid_full <- lipid_annot %>%
#   left_join(mm_baseline, by = "annotation")

# message("  Lipids with module assignment: ",
#         sum(!is.na(lipid_full$Module)), "/", nrow(lipid_full))

# # ── 4. Cross-validation report ────────────────────────────
# message("\nGenerating cross-validation report...")

# disagree_df <- lipid_full %>%
#   filter(CID_agree == "disagree") %>%
#   dplyr::select(annotation, InChiKey, Pubchem_CID_direct,
#                 Pubchem_CID_refmet, Super.class, Main.class)

# message("  CID disagreements: ", nrow(disagree_df))
# if (nrow(disagree_df) > 0) {
#   message("  These lipids have different PubChem CIDs in RefMet vs direct lookup.")
#   message("  The direct (InChiKey-based) CID is used as Pubchem_CID_final.")
#   message("  Review the validation file to manually verify if needed.")
# }

# write.csv(disagree_df,
#           file.path(out, "cid_disagreements.csv"),
#           row.names = FALSE)

# # ── 5. Class distribution by module ──────────────────────
# message("\nSummarizing lipid class distribution by module...")

# class_by_module <- lipid_full %>%
#   filter(!is.na(Module) & Module != "grey") %>%
#   group_by(Module, Main.class) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   arrange(Module, desc(n))

# write.csv(class_by_module,
#           file.path(out, "lipid_class_by_module.csv"),
#           row.names = FALSE)

# # ESI mode by module
# esi_by_module <- lipid_full %>%
#   filter(!is.na(Module) & Module != "grey") %>%
#   group_by(Module, ESI_mode) %>%
#   summarise(n = n(), .groups = "drop")

# write.csv(esi_by_module,
#           file.path(out, "esi_mode_by_module.csv"),
#           row.names = FALSE)

# # ── 6. Save full annotation table ─────────────────────────
# final_cols <- c("annotation", "Module", "identifier", "InChiKey",
#                 "ion_species", "mz", "ret_time", "ESI_mode",
#                 "Standardized name", "Formula", "Exact mass",
#                 "Super.class", "Main.class", "Sub.class",
#                 "Pubchem_CID_final", "Pubchem_CID_direct",
#                 "Pubchem_CID_refmet", "CID_agree",
#                 "CHEBI_ID", "HMDB_ID", "LM_ID", "KEGG_ID",
#                 "INCHI_KEY", "REFMET_ID")

# lipid_final <- lipid_full %>%
#   dplyr::select(any_of(final_cols)) %>%
#   arrange(Module, Super.class, Main.class, annotation)

# write.csv(lipid_final,
#           file.path(out, "lipid_annotations_full.csv"),
#           row.names = FALSE)
# write.xlsx(lipid_final,
#            file.path(out, "lipid_annotations_full.xlsx"))

# # Save as RDS for downstream use
# saveRDS(lipid_final,
#         file.path(cfg$output$processed,
#           paste0("lipid_annotations_", wgcna_run_tag, ".rds")))

# # ── 7. Hub lipid annotation tables ────────────────────────
# message("\nAnnotating hub lipids...")

# tp_labels <- c("Baseline", "Wk36_38", "Postpartum")
# for (tp in tp_labels) {
#   hub_path <- file.path("results/04_module_membership", wgcna_run_tag,
#                         paste0("hub_lipids_top10_", tp, ".csv"))
#   if (!file.exists(hub_path)) next

#   hubs <- read.csv(hub_path) %>%
#     rename(annotation = HubLipid) %>%
#     left_join(lipid_final %>%
#                 dplyr::select(annotation, Super.class, Main.class,
#                                Sub.class, Pubchem_CID_final,
#                                InChiKey, ESI_mode, mz, HMDB_ID),
#               by = "annotation")

#   write.xlsx(hubs,
#              file.path(out, paste0("hub_lipids_annotated_", tp, ".xlsx")))
#   message("  ", tp, ": ", nrow(hubs), " hub lipids annotated")
# }


# # ── 8. Create hub-based module ordering for heatmaps ──────
# # Goal: order modules biologically, with TG-associated modules first.
# # Uses annotated top hub lipids when available; falls back to all module lipids.

# message("\nCreating hub-based module order...")

# # Class group ordering — TGs always first, then DG, CE, PE/LPE, PC/LPC,
# # SM, Ceramides, Fatty acyls, Other. Edit the prefix numbers to reorder.
# class_group_from_text <- function(x) {
#   x <- ifelse(is.na(x), "", x)
#   x <- paste(x, collapse = " ")
#   dplyr::case_when(
#     grepl("Triradylglycerol|triacyl|triglycer|\\bTG\\b", x, ignore.case = TRUE) ~ "01_TG-associated",
#     grepl("Diradylglycerol|diacyl|diglycer|\\bDG\\b", x, ignore.case = TRUE) ~ "02_DG-associated",
#     grepl("cholest|\\bCE\\b", x, ignore.case = TRUE) ~ "03_Cholesteryl-esters",
#     grepl("Glycerophosphoethanolamine|phosphatidylethanolamine|\\bPE\\b|\\bLPE\\b", x, ignore.case = TRUE) ~ "04_PE/LPE",
#     grepl("Glycerophosphocholine|phosphatidylcholine|\\bPC\\b|\\bLPC\\b", x, ignore.case = TRUE) ~ "05_PC/LPC",
#     grepl("Phosphosphingolipid|sphingomyelin|\\bSM\\b", x, ignore.case = TRUE) ~ "06_Sphingomyelins",
#     grepl("Ceramide|\\bCer\\b", x, ignore.case = TRUE) ~ "07_Ceramides",
#     grepl("Fatty acid|Fatty ester|Fatty amide|fatty acyl|acylcarnitine", x, ignore.case = TRUE) ~ "08_Fatty-acyls",
#     grepl("unknown|other|unclassified", x, ignore.case = TRUE) ~ "99_Other/unknown",
#     TRUE ~ "90_Other"
#   )
# }

# hub_annot_all <- list()
# for (tp in tp_labels) {
#   hub_xlsx <- file.path(out, paste0("hub_lipids_annotated_", tp, ".xlsx"))
#   hub_csv  <- file.path(out, paste0("hub_lipids_annotated_", tp, ".csv"))
#   hub_raw  <- file.path("results/04_module_membership", wgcna_run_tag,
#                         paste0("hub_lipids_top10_", tp, ".csv"))

#   if (file.exists(hub_raw)) {
#     hubs_tp <- read.csv(hub_raw, stringsAsFactors = FALSE) %>%
#       rename(annotation = HubLipid) %>%
#       left_join(lipid_final %>%
#                   dplyr::select(annotation, Super.class, Main.class, Sub.class),
#                 by = "annotation") %>%
#       mutate(Timepoint = tp)
#     hub_annot_all[[tp]] <- hubs_tp
#     write.csv(hubs_tp, hub_csv, row.names = FALSE)
#   }
# }

# hub_annot_all <- bind_rows(hub_annot_all)

# if (nrow(hub_annot_all) > 0) {
#   module_class_summary <- hub_annot_all %>%
#     filter(!is.na(Module), Module != "grey") %>%
#     mutate(
#       ClassText = paste(Super.class, Main.class, Sub.class, sep = " | "),
#       ClassGroup = vapply(ClassText, class_group_from_text, character(1))
#     ) %>%
#     group_by(Module, ClassGroup) %>%
#     summarise(
#       n_hubs = n(),
#       mean_kME = mean(kME, na.rm = TRUE),
#       best_rank = min(Rank, na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     arrange(Module, desc(n_hubs), best_rank, desc(mean_kME)) %>%
#     group_by(Module) %>%
#     slice(1) %>%
#     ungroup()
# } else {
#   message("  No hub annotations found; falling back to all lipids per module.")
#   module_class_summary <- lipid_final %>%
#     filter(!is.na(Module), Module != "grey") %>%
#     mutate(
#       ClassText = paste(Super.class, Main.class, Sub.class, sep = " | "),
#       ClassGroup = vapply(ClassText, class_group_from_text, character(1))
#     ) %>%
#     group_by(Module, ClassGroup) %>%
#     summarise(n_hubs = n(), mean_kME = NA_real_, best_rank = NA_real_, .groups = "drop") %>%
#     arrange(Module, desc(n_hubs)) %>%
#     group_by(Module) %>%
#     slice(1) %>%
#     ungroup()
# }

# module_order_df <- module_class_summary %>%
#   mutate(
#     ClassOrder = as.integer(sub("_.*$", "", ClassGroup)),
#     Module = as.character(Module)
#   ) %>%
#   arrange(ClassOrder, ClassGroup, Module) %>%
#   mutate(PlotOrder = row_number()) %>%
#   dplyr::select(PlotOrder, Module, ClassGroup, n_hubs, mean_kME, best_rank)

# write.csv(module_order_df,
#           file.path(out, "module_order_by_hub_lipid_class.csv"),
#           row.names = FALSE)

# saveRDS(module_order_df,
#         file.path(cfg$output$processed,
#                   paste0("module_order_by_hub_lipid_class_", wgcna_run_tag, ".rds")))

# message("  Module order saved: module_order_by_hub_lipid_class.csv")
# print(module_order_df)


# # ── 9. Lipid class × module heatmaps (Figure 2) ──────────
# # heatmap_mode controls what is shown:
# #   "count"   → raw number of lipids per class per module
# #   "percent" → proportion of each module that is each class (%)
# #   "both"    → saves both versions
# heatmap_mode <- "both"   # ← change to "count" or "percent" if needed

# message("\nGenerating lipid class x module heatmaps (mode: ", heatmap_mode, ")...")

# # Build the count matrix: rows = lipid class, cols = module
# class_mod_mat <- lipid_full %>%
#   filter(!is.na(Module), Module != "grey", !is.na(Main.class)) %>%
#   group_by(Main.class, Module) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   tidyr::pivot_wider(names_from = Module, values_from = n, values_fill = 0)

# class_names <- class_mod_mat$Main.class
# count_mat   <- as.matrix(class_mod_mat[, -1])
# rownames(count_mat) <- class_names

# # Order modules by hub lipid class (biological ordering)
# if (file.exists(file.path(cfg$output$processed,
#     paste0("module_order_by_hub_lipid_class_", wgcna_run_tag, ".rds")))) {
#   mod_order <- readRDS(file.path(cfg$output$processed,
#     paste0("module_order_by_hub_lipid_class_", wgcna_run_tag, ".rds")))
#   ordered_mods <- mod_order$Module[mod_order$Module %in% colnames(count_mat)]
#   remaining    <- setdiff(colnames(count_mat), ordered_mods)
#   count_mat    <- count_mat[, c(ordered_mods, remaining), drop = FALSE]
# }

# # Order rows by total count descending
# row_totals <- rowSums(count_mat)
# count_mat  <- count_mat[order(row_totals, decreasing = TRUE), ]

# # Percent matrix
# pct_mat <- sweep(count_mat, 2, colSums(count_mat), "/") * 100

# # WGCNA module colors — used as column color bars at top of heatmap
# # R's built-in color names match WGCNA module names directly
# module_color_bar <- function(module_names) {
#   # Most WGCNA color names are valid R colors; handle exceptions
#   color_map <- setNames(module_names, module_names)
#   color_map["grey60"]       <- "grey60"
#   color_map["darkturquoise"]<- "darkturquoise"
#   color_map["darkolivegreen"]<- "darkolivegreen"
#   color_map["saddlebrown"]  <- "saddlebrown"
#   color_map["paleturquoise"]<- "paleturquoise"
#   color_map["yellowgreen"]  <- "yellowgreen"
#   color_map["steelblue"]    <- "steelblue"
#   color_map["darkmagenta"]  <- "darkmagenta"
#   color_map["darkorange"]   <- "darkorange"
#   color_map["darkred"]      <- "darkred"
#   color_map["darkgreen"]    <- "darkgreen"
#   color_map["darkgrey"]     <- "darkgrey"
#   color_map["sienna3"]      <- "sienna3"
#   color_map["violet"]       <- "violet"
#   color_map["white"]        <- "grey90"   # white module → light grey so visible
#   color_map
# }

# # ── Display name mapping (figures only — data tables unchanged) ───
# # Maps RefMet Main.class names to reader-friendly labels for plots
# # Excel/CSV outputs always use original RefMet names
# LIPID_CLASS_DISPLAY <- c(
#   "Triradylglycerols"                        = "Triglycerides (TG)",
#   "Diradylglycerols"                         = "Diglycerides (DG)",
#   "Monoradylglycerols"                        = "Monoglycerides (MG)",
#   "Glycerophosphocholines"                   = "Phosphatidylcholines (PC/LPC)",
#   "Glycerophosphoethanolamines"              = "Phosphatidylethanolamines (PE/LPE)",
#   "Glycerophosphoinositols"                  = "Phosphatidylinositols (PI/LPI)",
#   "Glycerophosphoserines"                    = "Phosphatidylserines (PS)",
#   "Glycerophosphatidylethanol"               = "Phosphatidylethanol (PEtOH)",
#   "Phosphosphingolipids"                     = "Sphingomyelins (SM)",
#   "Neutral_glycosphingolipids"               = "Neutral Glycosphingolipids",
#   "Glycosyldiradylglycerols"                 = "Glycosyldiacylglycerols",
#   "Perfluoroalkyl_Phosphate_Acid_derivative" = "PFAS Derivatives",
#   "Ceramides"                                = "Ceramides",
#   "Glycosphingolipids"                       = "Glycosphingolipids",
#   "Sphingoid bases"                          = "Sphingoid Bases",
#   "Fatty acids"                              = "Fatty Acids",
#   "Fatty esters"                             = "Fatty Esters",
#   "Fatty amides"                             = "Fatty Amides",
#   "Eicosanoids"                              = "Eicosanoids",
#   "Docosanoids"                              = "Docosanoids",
#   "Octadecanoids"                            = "Octadecanoids",
#   "Sterol esters"                            = "Sterol Esters",
#   "Sterols"                                  = "Sterols",
#   "Glycerophosphates"                        = "Glycerophosphates",
#   "Keto acids"                               = "Keto Acids"
# )

# # Helper: apply display names to a vector of class names
# apply_display_names <- function(x) {
#   ifelse(x %in% names(LIPID_CLASS_DISPLAY),
#          LIPID_CLASS_DISPLAY[x],
#          x)   # fallback: keep original if not in map
# }

# make_class_heatmap <- function(mat, title_str, legend_label,
#                                 fmt_fn = function(x) ifelse(x == 0, "", as.character(x))) {
#   max_val    <- max(mat, na.rm = TRUE)
#   n_rows     <- nrow(mat)
#   n_cols     <- ncol(mat)
#   mod_names  <- colnames(mat)
#   mod_colors <- module_color_bar(mod_names)

#   # ── Main tile data ────────────────────────────────────────
#   # Apply display names to row labels (figures only)
#   display_rownames <- apply_display_names(rownames(mat))
#   mat_display <- mat
#   rownames(mat_display) <- display_rownames

#   df <- as.data.frame(mat_display) %>%
#     tibble::rownames_to_column("LipidClass") %>%
#     tidyr::pivot_longer(-LipidClass, names_to = "Module", values_to = "Value") %>%
#     mutate(
#       LipidClass  = factor(LipidClass, levels = rev(display_rownames)),
#       Module      = factor(Module, levels = mod_names),
#       Label       = fmt_fn(Value),
#       TextColor   = ifelse(Value > max_val * 0.55, "white", "black"),
#       # Zebra stripe: alternating row backgrounds
#       RowIdx      = as.integer(LipidClass),
#       RowStripe   = ifelse(RowIdx %% 2 == 0, "#F7F7F7", "#FFFFFF")
#     )

#   # ── Module color bar data (row above x-axis labels) ───────
#   color_bar_df <- data.frame(
#     Module    = factor(mod_names, levels = mod_names),
#     FillColor = unname(mod_colors[mod_names]),
#     y         = "Module color"
#   )

#   # ── Main plot ─────────────────────────────────────────────
#   p <- ggplot(df, aes(x = Module, y = LipidClass)) +
#     # Zebra stripe background (full row width)
#     geom_tile(aes(fill = RowStripe), alpha = 0.6, colour = NA) +
#     scale_fill_identity() +
#     ggnewscale::new_scale_fill() +
#     # Value tiles
#     geom_tile(aes(fill = Value), colour = NA,
#               width = 0.95, height = 0.95) +
#     # Vertical grid lines between every module column
#     geom_vline(
#       xintercept = seq(0.5, n_cols + 0.5, by = 1),
#       colour = "grey80", linewidth = 0.25
#     ) +
#     # Horizontal grid lines between every lipid class row
#     geom_hline(
#       yintercept = seq(0.5, n_rows + 0.5, by = 1),
#       colour = "grey80", linewidth = 0.25
#     ) +
#     # Cell text
#     geom_text(aes(label = Label, colour = TextColor),
#               size = FS_SMALL / ggplot2::.pt, fontface = "bold") +
#     scale_colour_identity() +
#     scale_fill_gradientn(
#       colours  = c("#FFFFFF", "#FDE8D8", "#F4A261", "#E63946", "#8B0000"),
#       values   = scales::rescale(c(0, 0.2, 0.5, 0.75, 1)),
#       na.value = "grey97",
#       name     = legend_label
#     ) +
#     # Outer border
#     annotate("rect",
#              xmin = 0.5, xmax = n_cols + 0.5,
#              ymin = 0.5, ymax = n_rows + 0.5,
#              fill = NA, colour = "grey50", linewidth = 0.4) +
#     labs(title = title_str, x = NULL, y = "Lipid class") +
#     scale_x_discrete(expand = c(0, 0)) +
#     scale_y_discrete(expand = c(0, 0)) +
#     theme_growell() +
#     theme(
#       axis.text.x     = element_text(angle = 50, hjust = 1,
#                                      size  = FS_SMALL - 0.5, colour = "black"),
#       axis.text.y     = element_text(size = FS_SMALL, colour = "black"),
#       axis.title.y    = element_text(size = FS_BASE),
#       axis.ticks      = element_blank(),
#       legend.position = "right",
#       panel.border    = element_blank(),
#       panel.grid      = element_blank()
#     )

#   # ── Module color bar aligned to heatmap columns ───────────
#   p_bar <- ggplot(color_bar_df, aes(x = Module, y = y, fill = FillColor)) +
#     geom_tile(colour = "white", linewidth = 0.5, height = 0.85) +
#     scale_fill_identity() +
#     scale_x_discrete(limits = mod_names, expand = c(0, 0)) +
#     labs(x = "Module", y = NULL) +
#     coord_cartesian(clip = "off") +
#     theme_void() +
#     theme(
#       axis.text.x  = element_blank(),
#       axis.title.x = element_text(size = FS_BASE, vjust = -0.5),
#       plot.margin  = unit(c(2, 0, 0, 0), "pt")
#     )

#   # Stack: color bar on top of heatmap, aligned
#   p_bar / p +
#     patchwork::plot_layout(heights = c(0.05, 1)) &
#     theme(plot.margin = unit(c(1, 1, 1, 1), "pt"))
# }

# # Save based on mode
# save_heatmap <- function(p, suffix, panel_id) {
#   fname <- file.path(out, paste0("heatmap_lipid_class_module_", suffix, ".pdf"))
#   ggsave(fname, plot = p,
#          width  = FIG_WIDTH_FULL, height = 120, units = "mm",
#          device = cairo_pdf)
#   ggsave(gsub(".pdf", ".png", fname), plot = p,
#          width = FIG_WIDTH_FULL, height = 120, units = "mm",
#          dpi   = FIGURE_DPI)
#   save_panel(p, "fig2", panel_id,
#              width_mm = FIG_WIDTH_FULL, height_mm = 120)
#   message("  Saved: heatmap_lipid_class_module_", suffix, ".pdf → fig2_", panel_id)
# }

# if (heatmap_mode %in% c("count", "both")) {
#   p_count <- make_class_heatmap(
#     count_mat,
#     title_str    = "Lipid class distribution by module (count)",
#     legend_label = "No. lipids",
#     fmt_fn       = function(x) ifelse(x == 0, "", as.character(x))
#   )
#   save_heatmap(p_count, "count", if (heatmap_mode == "both") "D" else "D")
# }

# if (heatmap_mode %in% c("percent", "both")) {
#   p_pct <- make_class_heatmap(
#     pct_mat,
#     title_str    = "Lipid class distribution by module (%)",
#     legend_label = "% of module",
#     fmt_fn       = function(x) ifelse(x < 1, "", paste0(round(x), "%"))
#   )
#   save_heatmap(p_pct, "percent", if (heatmap_mode == "both") "E" else "D")
# }

# # ── 10. Annotated lipid overview bar chart (Figure 2A) ────
# # Bar chart showing number of identified lipids per Main class
# # coloured by Super class — overview of the full lipidome detected
# message("\nGenerating annotated lipid overview chart...")

# # Colorblind-safe palette for lipid super classes
# superclass_colors <- c(
#   "Fatty Acyls"            = "#E63946",
#   "Glycerophospholipids"   = "#1D3D8F",
#   "Sphingolipids"          = "#F4A261",
#   "Sterol Lipids"          = "#2A9D8F",
#   "Glycerolipids"          = "#E9C46A",
#   "Prenol Lipids"          = "#264653",
#   "Saccharolipids"         = "#A8DADC",
#   "Polyketides"            = "#457B9D",
#   "Other"                  = "#CCCCCC"
# )

# lipid_overview <- lipid_full %>%
#   mutate(
#     Super.class = ifelse(is.na(Super.class) | Super.class == "",
#                           "Other", Super.class),
#     Main.class  = ifelse(is.na(Main.class) | Main.class == "",
#                           "Unknown", Main.class)
#   ) %>%
#   group_by(Super.class, Main.class) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   arrange(Super.class, desc(n)) %>%
#   mutate(
#     # Apply display names before making factor so levels are correct
#     Main.class_display = apply_display_names(as.character(Main.class)),
#     Main.class_display = factor(Main.class_display,
#                                  levels = rev(unique(Main.class_display)))
#   )

# p_overview <- ggplot(lipid_overview,
#                       aes(x = n, y = Main.class_display, fill = Super.class)) +
#   geom_col(width = 0.7) +
#   scale_fill_manual(values = superclass_colors,
#                     name   = "Lipid super class") +
#   labs(title = "Identified lipids by class",
#        x = "Number of lipids", y = NULL) +
#   theme_growell(grid = "y") +
#   theme(legend.position = "right",
#         legend.text     = element_text(size = FS_SMALL),
#         axis.text.y     = element_text(size = FS_SMALL))

# ggsave(file.path(out, "annotated_lipid_overview.pdf"),
#        plot = p_overview,
#        width = FIG_WIDTH_FULL, height = 100, units = "mm",
#        device = cairo_pdf)

# save_panel(p_overview, "fig2", "A",
#            width_mm = FIG_WIDTH_FULL, height_mm = 100)
# message("  Saved: annotated_lipid_overview.pdf → fig2_A")

# message("\nScript 07 complete → ", out)
# message("  Full annotation: lipid_annotations_full.xlsx")
# message("  CID disagreements to review: cid_disagreements.csv (", nrow(disagree_df), " lipids)")