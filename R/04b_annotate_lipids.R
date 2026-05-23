# R/04b_annotate_lipids.R
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
# Outputs → results/04b_annotations/<run_tag>/
# -------------------------------------------------------

source("R/00_load_packages.R")
cfg <- yaml::read_yaml("config/config.yml")

wgcna_run_tag <- trimws(readLines(
  file.path(cfg$output$processed, "current_wgcna_run_tag.txt")))
message("WGCNA run: ", wgcna_run_tag)

out <- file.path("results/04b_annotations", wgcna_run_tag)
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

class_group_from_text <- function(x) {
  x <- ifelse(is.na(x), "", x)
  x <- paste(x, collapse = " ")
  dplyr::case_when(
    grepl("triacyl|triglycer|\\bTG\\b", x, ignore.case = TRUE) ~ "01_TG-associated",
    grepl("diacyl|diglycer|\\bDG\\b", x, ignore.case = TRUE) ~ "02_DG-associated",
    grepl("cholest|\\bCE\\b", x, ignore.case = TRUE) ~ "03_Cholesteryl-esters",
    grepl("phosphatidylcholine|glycerophosphocholine|\\bPC\\b|lysophosphatidylcholine|\\bLPC\\b", x, ignore.case = TRUE) ~ "04_PC/LPC",
    grepl("phosphatidylethanolamine|\\bPE\\b|lysophosphatidylethanolamine|\\bLPE\\b", x, ignore.case = TRUE) ~ "05_PE/LPE",
    grepl("sphingomyelin|\\bSM\\b", x, ignore.case = TRUE) ~ "06_Sphingomyelins",
    grepl("ceramide|\\bCer\\b", x, ignore.case = TRUE) ~ "07_Ceramides",
    grepl("fatty acyl|fatty acid|acylcarnitine", x, ignore.case = TRUE) ~ "08_Fatty-acyls",
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

message("\nScript 04b complete → ", out)
message("  Full annotation: lipid_annotations_full.xlsx")
message("  CID disagreements to review: cid_disagreements.csv (", nrow(disagree_df), " lipids)")
# # R/04b_annotate_lipids.R
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
# # Outputs → results/04b_annotations/<run_tag>/
# # -------------------------------------------------------

# source("R/00_load_packages.R")
# cfg <- yaml::read_yaml("config/config.yml")

# wgcna_run_tag <- trimws(readLines(
#   file.path(cfg$output$processed, "current_wgcna_run_tag.txt")))
# message("WGCNA run: ", wgcna_run_tag)

# out <- file.path("results/04b_annotations", wgcna_run_tag)
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

# message("\nScript 04b complete → ", out)
# message("  Full annotation: lipid_annotations_full.xlsx")
# message("  CID disagreements to review: cid_disagreements.csv (", nrow(disagree_df), " lipids)")