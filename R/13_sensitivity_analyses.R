# R/13_sensitivity_analyses.R
# -------------------------------------------------------
# Sensitivity analyses (reviewer request)
# Outputs → results/10_sensitivity/
# -------------------------------------------------------
source("R/00_load_packages.R")
cfg <- yaml::read_yaml("config/config.yml")

# Read analysis tags set by scripts 03 and 05
wgcna_run_tag <- trimws(readLines(
  file.path(cfg$output$processed, "current_wgcna_run_tag.txt")))
analysis_tag_raw <- trimws(readLines(
  file.path(cfg$output$processed, "current_analysis_tag.txt")))
analysis_parts <- strsplit(analysis_tag_raw, "\\|")[[1]]
wgcna_run_tag  <- analysis_parts[1]
method         <- analysis_parts[2]
rds_tag        <- paste0(wgcna_run_tag, "_", method)
message("WGCNA run: ", wgcna_run_tag, " | Method: ", method)

out <- file.path(cfg$output$s13, wgcna_run_tag, method)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

consensusMEs_ID      <- readRDS(file.path(cfg$output$processed, paste0("consensusMEs_ID_", wgcna_run_tag, ".rds")))
moduleTraitCor_ID    <- readRDS(file.path(cfg$output$processed, paste0("moduleTraitCor_ID_", rds_tag, ".rds")))
moduleTraitPvalue_ID <- readRDS(file.path(cfg$output$processed, paste0("moduleTraitPvalue_ID_", rds_tag, ".rds")))
demographic_data     <- readRDS(file.path(cfg$output$processed, "demographic_data_bmi.rds"))
exprSize_ID          <- readRDS(file.path(cfg$output$processed, "exprSize_ID.rds"))
All_nSets_ID <- exprSize_ID$nSets
tp_labels    <- c("Baseline","Wk36_38","Postpartum")
key_traits   <- c("apo","ppwr_e","apo_hdp","apo_gdm","bmi_base","gwg")
key_modules  <- c("magenta","purple","brown","tan")

cor_pval <- function(me_df, tr_df) {
  ids <- intersect(rownames(me_df), rownames(tr_df))
  me  <- me_df[ids,]; tr <- tr_df[ids,]
  cor <- WGCNA::cor(me, tr, method=if(method=="spearman") "spearman" else "pearson", use="p")
  pv  <- WGCNA::corPvalueFisher(cor, nrow(me), twoSided=TRUE)
  list(cor=cor, pval=pv, n=nrow(me))
}
save_wb <- function(res, path) {
  wb <- createWorkbook()
  addWorksheet(wb,"Correlations"); addWorksheet(wb,"P_values")
  writeData(wb,"Correlations", round(res$cor,3),  rowNames=TRUE)
  writeData(wb,"P_values",     signif(res$pval,3), rowNames=TRUE)
  saveWorkbook(wb, path, overwrite=TRUE)
}

# ── 1. Complete-case ─────────────────────────────────────
message("Sensitivity 1: Complete-case analysis...")
for (set in 1:All_nSets_ID) {
  tp  <- tp_labels[set]
  me  <- consensusMEs_ID[[set]]$data
  tr  <- demographic_data[, intersect(key_traits, colnames(demographic_data))]
  ids <- intersect(rownames(me), rownames(tr))
  cc  <- complete.cases(cbind(me[ids,], tr[ids,]))
  message(sprintf("  %s: %d complete-case samples", tp, sum(cc)))
  res <- cor_pval(me[ids[cc],], tr[ids[cc],])
  save_wb(res, file.path(out, paste0("sens1_complete_case_",tp,".xlsx")))
}

# ── 2. Intervention arm stratification ───────────────────
message("Sensitivity 2: Arm stratification...")
if ("group" %in% colnames(demographic_data)) {
  for (arm in c(0,1)) {
    lbl <- if(arm==0) "Control" else "Intervention"
    ids_arm <- rownames(demographic_data)[demographic_data$group==arm]
    for (set in 1:All_nSets_ID) {
      tp  <- tp_labels[set]; me <- consensusMEs_ID[[set]]$data
      ids <- intersect(ids_arm, rownames(me)); if (length(ids)<5) next
      tr  <- demographic_data[ids, intersect(key_traits, colnames(demographic_data))]
      res <- cor_pval(me[ids,], tr)
      message(sprintf("  %s — %s: n=%d", lbl, tp, res$n))
      save_wb(res, file.path(out, paste0("sens2_",lbl,"_",tp,".xlsx")))
    }
  }
}

# ── 3. Outlier removal ───────────────────────────────────
message("Sensitivity 3: Outlier removal (>3 SD on PC1)...")
for (set in 1:All_nSets_ID) {
  tp  <- tp_labels[set]; me <- consensusMEs_ID[[set]]$data
  pc1 <- prcomp(me, scale.=TRUE)$x[,1]
  out_ids <- names(pc1[abs(pc1-mean(pc1)) > 3*sd(pc1)])
  message(sprintf("  %s: removed %d outliers", tp, length(out_ids)))
  me_clean <- me[!rownames(me) %in% out_ids,]
  tr_clean <- demographic_data[rownames(me_clean),
               intersect(key_traits, colnames(demographic_data))]
  res <- cor_pval(me_clean, tr_clean)
  save_wb(res, file.path(out, paste0("sens3_no_outliers_",tp,".xlsx")))
}

# ── 4. BMI sensitivity ───────────────────────────────────
message("Sensitivity 4: BMI as confounder vs mediator...")
non_bmi <- setdiff(key_traits, c("bmi_base","bmi"))
for (set in 1:All_nSets_ID) {
  tp  <- tp_labels[set]; me <- consensusMEs_ID[[set]]$data
  ids <- intersect(rownames(me), rownames(demographic_data))
  rows <- lapply(colnames(me), function(m) lapply(non_bmi, function(tr) {
    if (!all(c(tr,"bmi_base") %in% colnames(demographic_data))) return(NULL)
    df <- na.omit(data.frame(ME=me[ids,m], Trait=demographic_data[ids,tr],
                              BMI=demographic_data[ids,"bmi_base"]))
    if (nrow(df)<10) return(NULL)
    u <- summary(lm(ME~Trait,     data=df))$coefficients
    a <- summary(lm(ME~Trait+BMI, data=df))$coefficients
    data.frame(ME=m, Trait=tr, Timepoint=tp,
               beta_unadj=u["Trait","Estimate"], p_unadj=u["Trait","Pr(>|t|)"],
               beta_adj  =a["Trait","Estimate"], p_adj  =a["Trait","Pr(>|t|)"],
               n=nrow(df))
  }) %>% bind_rows()) %>% bind_rows()
  write.xlsx(rows, file.path(out, paste0("sens4_bmi_adjustment_",tp,".xlsx")))
}

# ── 5. Cross-timepoint consistency ───────────────────────
message("Sensitivity 5: Cross-timepoint consistency...")
consistency <- lapply(paste0("ME",key_modules), function(me)
  lapply(key_traits, function(tr) {
    rs <- sapply(1:All_nSets_ID, function(s) {
      m <- moduleTraitCor_ID[[s]]
      if (me %in% rownames(m) && tr %in% colnames(m)) m[me,tr] else NA })
    ps <- sapply(1:All_nSets_ID, function(s) {
      m <- moduleTraitPvalue_ID[[s]]
      if (me %in% rownames(m) && tr %in% colnames(m)) m[me,tr] else NA })
    data.frame(Module=gsub("ME","",me), Trait=tr,
               r_Base=rs[1], p_Base=ps[1], r_Wk36=rs[2], p_Wk36=ps[2],
               r_PP=rs[3], p_PP=ps[3],
               consistent=sum(!is.na(ps)&ps<0.05)>=2)
  }) %>% bind_rows()) %>% bind_rows()
write.xlsx(consistency, file.path(out, "sens5_cross_timepoint.xlsx"))

pdf(file.path(out, "sens5_cross_timepoint.pdf"), width=10, height=6)
consistency %>%
  tidyr::pivot_longer(cols=c(r_Base,r_Wk36,r_PP), names_to="Timepoint", values_to="r") %>%
  mutate(Timepoint=gsub("r_","",Timepoint), label=paste(Module,Trait,sep=" x ")) %>%
  ggplot(aes(x=Timepoint, y=r, group=label, color=Module)) +
  geom_line(alpha=0.7) + geom_point(size=2) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  facet_wrap(~Trait, scales="free_y") +
  labs(title="Cross-timepoint consistency", y="Pearson r") +
  theme_bw(base_size=11) + theme(panel.grid=element_blank()) %>% print()
dev.off()

message("Script 13 complete → ", out)