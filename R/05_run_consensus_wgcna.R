# R/05_run_consensus_wgcna.R
# -------------------------------------------------------
# Soft power selection and consensus WLCNA
# Outputs → results/03_consensus_wgcna/
# -------------------------------------------------------
source("R/00_load_packages.R")
cfg <- yaml::read_yaml("config/config.yml")

wgcna_run_tag <- sprintf("power%d_split%d_merge%03d",
  cfg$wgcna$soft_power,
  cfg$wgcna$deep_split,
  as.integer(cfg$wgcna$merge_cut_height * 100))
message("WGCNA run tag: ", wgcna_run_tag)

out <- file.path(cfg$output$s05, wgcna_run_tag)
dir.create(out, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$output$processed, showWarnings = FALSE, recursive = TRUE)

multiExpr_ID <- readRDS(file.path(cfg$output$processed, "multiExpr_ID.rds"))
exprSize_ID  <- readRDS(file.path(cfg$output$processed, "exprSize_ID.rds"))
nSets_ID     <- exprSize_ID$nSets
setLabels_ID <- cfg$wgcna$timepoint_labels

message("Running soft-threshold analysis...")
powers_ID      <- c(seq(4,10,by=1), seq(12,26,by=2))
powerTables_ID <- vector(mode="list", length=nSets_ID)
for (set in 1:nSets_ID) {
  powerTables_ID[[set]] <- list(data = pickSoftThreshold(
    multiExpr_ID[[set]]$data, blockSize=1000, powerVector=powers_ID,
    corFnc=cfg$wgcna$cor_type, networkType=cfg$wgcna$network_type, verbose=2)[[2]])
}
collectGarbage()

colors_ID <- c("black","red","green")
plotCols_ID <- c(2,5,6,7)
colNames_ID <- c("Scale Free Topology Model Fit","Mean connectivity",
                 "Median connectivity","Max connectivity")
ylim_ID <- matrix(NA, nrow=2, ncol=4)
for (set in 1:nSets_ID) for (col in 1:4) {
  ylim_ID[1,col] <- min(ylim_ID[1,col], powerTables_ID[[set]]$data[,plotCols_ID[col]], na.rm=TRUE)
  ylim_ID[2,col] <- max(ylim_ID[2,col], powerTables_ID[[set]]$data[,plotCols_ID[col]], na.rm=TRUE)
}
pdf(file.path(out, "scaleFreeAnalysis.pdf"), width=8, height=6)
par(mfcol=c(2,2)); par(mar=c(4.2,4.2,2.2,0.5)); cex1 <- 0.7
for (col in 1:4) for (set in 1:nSets_ID) {
  if (set==1) {
    plot(powerTables_ID[[set]]$data[,1],
         -sign(powerTables_ID[[set]]$data[,3])*powerTables_ID[[set]]$data[,2],
         xlab="Soft Threshold (power)", ylab=colNames_ID[col],
         type="n", ylim=ylim_ID[,col], main=colNames_ID[col])
    abline(h=c(0.80,0.90), col="red"); addGrid()
    abline(v=cfg$wgcna$soft_power, col="blue", lty=2, lwd=1.5)
    text(cfg$wgcna$soft_power, ylim_ID[2,col]*0.95,
         paste0("power=", cfg$wgcna$soft_power),
         col="blue", cex=0.7, adj=c(-0.1,1))
  }
  if (col==1) {
    text(powerTables_ID[[set]]$data[,1],
         -sign(powerTables_ID[[set]]$data[,3])*powerTables_ID[[set]]$data[,2],
         labels=powers_ID, cex=cex1, col=colors_ID[set])
  } else {
    text(powerTables_ID[[set]]$data[,1], powerTables_ID[[set]]$data[,plotCols_ID[col]],
         labels=powers_ID, cex=cex1, col=colors_ID[set])
  }
  legend(if(col==1)"bottomright" else "topright",
         legend=setLabels_ID, col=colors_ID, pch=20)
}
dev.off()

message("Running blockwiseConsensusModules...")
consensusMods_ID <- blockwiseConsensusModules(
  multiExpr_ID, checkMissingData=FALSE,
  maxBlockSize       = cfg$wgcna$max_block_size,
  corType            = cfg$wgcna$cor_type,
  maxPOutliers       = cfg$wgcna$max_p_outliers,
  power              = cfg$wgcna$soft_power,
  networkType        = cfg$wgcna$network_type,
  checkPower         = FALSE,
  TOMType            = cfg$wgcna$network_type,
  networkCalibration = "full quantile",
  saveConsensusTOMs  = TRUE,
  consensusTOMFilePattern = file.path(cfg$output$processed,
    paste0("consensusTOM_", wgcna_run_tag, "-block.%b.RData")),
  deepSplit          = cfg$wgcna$deep_split,
  minModuleSize      = cfg$wgcna$min_module_size,
  mergeCutHeight     = cfg$wgcna$merge_cut_height,
  verbose            = 5)

module_dist <- as.data.frame(table(consensusMods_ID$colors) %>% sort(decreasing=TRUE))
colnames(module_dist) <- c("Module","Lipids")
write.csv(module_dist, file.path(out, "module_distribution.csv"), row.names=FALSE)
message("Modules: ", nrow(module_dist)); print(module_dist)

pdf(file.path(out, "cluster_dendrogram.pdf"), width=10, height=5)
plotDendroAndColors(dendro=consensusMods_ID$dendrograms[[1]],
  colors=consensusMods_ID$colors, groupLabels="Modules",
  dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05,
  marAll=c(1,5,1,0), main=wgcna_run_tag, cex.colorLabels=1.3)
dev.off()

tp_names <- c("Identified_Baseline","Identified_TP36_38weeks","Identified_Postpartum")
tp_short <- c("Baseline","Wk36_38","Postpartum")
for (i in seq_along(tp_names)) {
  ME_data <- consensusMods_ID$multiMEs[[tp_names[i]]]$data
  METree  <- (1 - bicor(ME_data, maxPOutliers=0.1)) %>% as.dist() %>% hclust(method="average")
  pdf(file.path(out, paste0("eigennode_dendrogram_",tp_short[i],".pdf")), height=5, width=10)
  par(mar=c(0,5,1,1))
  plot(METree, main="", xlab="", sub="", ylim=c(0,1), cex=0.6)
  abline(h=cfg$wgcna$merge_cut_height, col="red"); dev.off()
}

consensusMEs_ID <- consensusOrderMEs(consensusMods_ID$multiMEs)
pdf(file.path(out, "eigennode_network.pdf"), width=8, height=7)
par(cex=0.8)
plotEigengeneNetworks(consensusMEs_ID, setLabels=setLabels_ID,
  plotDendrograms=FALSE, marHeatmap=c(3,3,2,1),
  zlimPreservation=c(0.5,1), xLabelsAngle=90)
dev.off()

saveRDS(consensusMods_ID, file.path(cfg$output$processed,
  paste0("consensusMods_ID_", wgcna_run_tag, ".rds")))
saveRDS(consensusMEs_ID, file.path(cfg$output$processed,
  paste0("consensusMEs_ID_", wgcna_run_tag, ".rds")))
writeLines(wgcna_run_tag,
  file.path(cfg$output$processed, "current_wgcna_run_tag.txt"))

message("Script 05 complete → ", out)
