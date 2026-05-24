# R/12_plsda.R
# -------------------------------------------------------
# PLS-DA with permutation testing (reviewer request)
# Outputs → results/09_plsda/
# -------------------------------------------------------
source("R/00_load_packages.R")
source("R/00_figure_theme.R")
cfg <- yaml::read_yaml("config/config.yml")
out <- cfg$output$s12
dir.create(out, showWarnings = FALSE, recursive = TRUE)

multiExpr_ID  <- readRDS(file.path(cfg$output$processed, "multiExpr_ID.rds"))
status_labels <- readRDS(file.path(cfg$output$processed, "lipid_status_labels.rds"))
tp_names  <- names(multiExpr_ID); tp_labels <- c("Baseline","Wk36_38","Postpartum")
outcomes  <- c("apo_status","ppwr_status","egwg_status")
ncomp     <- cfg$plsda$ncomp; folds <- cfg$plsda$folds; nrepeat <- cfg$plsda$nrepeat

run_plsda <- function(expr, group_vec, outcome_name, tp_label) {
  common <- intersect(rownames(expr), names(group_vec))
  X <- expr[common,]; Y <- factor(group_vec[common])
  valid_cols <- apply(X,2,function(x) !all(is.na(x)) && var(x,na.rm=TRUE)>0)
  X <- X[,valid_cols]; cc <- complete.cases(X,Y); X <- X[cc,]; Y <- Y[cc]
  if (length(unique(Y))<2 || nrow(X)<10) { message("  Skipping: insufficient data"); return(NULL) }
  message(sprintf("  PLS-DA: %s x %s (n=%d, p=%d)", tp_label, outcome_name, nrow(X), ncol(X)))
  fit <- tryCatch(mixOmics::plsda(X, Y, ncomp=ncomp), error=function(e) { message(e$message); NULL })
  if (is.null(fit)) return(NULL)
  message("  Permutation testing (", nrepeat, " repeats)...")
  perf <- tryCatch(
    mixOmics::perf(fit, validation=cfg$plsda$validation, folds=folds, nrepeat=nrepeat, progressBar=FALSE),
    error=function(e) { message("  Permutation error: ", e$message); NULL })
  pfx <- file.path(out, paste0("plsda_", tp_label, "_", outcome_name))
  pdf(paste0(pfx, "_scoreplot.pdf"), width=7, height=6)
  mixOmics::plotIndiv(fit, comp=c(1,2), legend=TRUE, ind.names=FALSE, ellipse=TRUE,
    title=paste("PLS-DA:", tp_label, "\u2014", outcome_name)); dev.off()
  pdf(paste0(pfx, "_loadings.pdf"), width=8, height=6)
  mixOmics::plotLoadings(fit, comp=1, contrib="max", method="mean", ndisplay=20,
    title=paste("Loadings Comp1:", tp_label, "\u2014", outcome_name)); dev.off()
  if (!is.null(perf)) {
    pdf(paste0(pfx, "_performance.pdf"), width=7, height=5)
    plot(perf, col=color.mixo(1:3), sd=TRUE, legend.position="horizontal"); dev.off()
  }
  vip_df <- data.frame(Lipid=rownames(mixOmics::vip(fit)),
                        VIP_comp1=mixOmics::vip(fit)[,1]) %>% arrange(desc(VIP_comp1))
  wb <- createWorkbook(); addWorksheet(wb,"VIP_scores"); writeData(wb,"VIP_scores", vip_df)
  if (!is.null(perf)) {
    addWorksheet(wb,"Performance")
    writeData(wb,"Performance", as.data.frame(perf$error.rate$overall))
  }
  saveWorkbook(wb, paste0(pfx, ".xlsx"), overwrite=TRUE)
  list(plsda=fit, perf=perf, vip=vip_df, tp=tp_label, outcome=outcome_name,
       n_samples=nrow(X), n_lipids=ncol(X))
}

plsda_results <- list()
for (i in seq_along(tp_names)) {
  expr <- multiExpr_ID[[tp_names[i]]]$data; tp <- tp_labels[i]
  for (outcome in outcomes) {
    if (!(outcome %in% colnames(status_labels))) next
    grp_vec <- setNames(status_labels[[outcome]], rownames(status_labels))
    plsda_results[[paste0(tp,"_",outcome)]] <- run_plsda(expr, grp_vec, outcome, tp)
  }
}
saveRDS(plsda_results, file.path(cfg$output$processed, "plsda_results.rds"))

summary_rows <- lapply(names(plsda_results), function(key) {
  res <- plsda_results[[key]]
  if (is.null(res)) return(data.frame(Key=key, Status="Failed"))
  data.frame(Key=key, Status="Complete", N_samples=res$n_samples,
             N_lipids=res$n_lipids, Top_VIP=res$vip$Lipid[1],
             VIP_val=round(res$vip$VIP_comp1[1],3))
}) %>% bind_rows()
write.xlsx(summary_rows, file.path(out, "plsda_summary.xlsx"))
print(summary_rows)
message("Script 12 complete → ", out)