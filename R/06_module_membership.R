# R/06_module_membership.R
# -------------------------------------------------------
# kME computation and hub lipid identification
# Outputs → results/04_module_membership/<run_tag>/
# -------------------------------------------------------
source("R/00_load_packages.R")
cfg <- yaml::read_yaml("config/config.yml")

wgcna_run_tag <- trimws(readLines(
  file.path(cfg$output$processed, "current_wgcna_run_tag.txt")))
message("WGCNA run: ", wgcna_run_tag)

out <- file.path(cfg$output$s06, wgcna_run_tag)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

multiExpr_ID     <- readRDS(file.path(cfg$output$processed, "multiExpr_ID.rds"))
consensusMods_ID <- readRDS(file.path(cfg$output$processed,
  paste0("consensusMods_ID_", wgcna_run_tag, ".rds")))
consensusMEs_ID  <- readRDS(file.path(cfg$output$processed,
  paste0("consensusMEs_ID_", wgcna_run_tag, ".rds")))

tp_names  <- c("Identified_Baseline","Identified_TP36_38weeks","Identified_Postpartum")
tp_labels <- c("Baseline","Wk36_38","Postpartum")

message("Computing module membership...")
moduleMembership_ID <- mtd.mapply(bicorAndPvalue, multiExpr_ID, consensusMEs_ID,
  MoreArgs = list(alternative="two.sided", use="pairwise.complete.obs",
                  maxPOutliers=cfg$wgcna$max_p_outliers))

mm_list <- list()
for (i in seq_along(tp_names)) {
  mm <- as.data.frame(moduleMembership_ID[[tp_names[i]]]$data$bicor)
  colnames(mm) <- gsub("ME","",colnames(mm))
  mm$Probe  <- rownames(mm)
  mm$Module <- consensusMods_ID$colors
  mm_list[[tp_labels[i]]] <- mm
  write.xlsx(mm, file.path(out, paste0("module_membership_", tp_labels[i], ".xlsx")))
  write.table(mm, file.path(out, paste0("module_membership_", tp_labels[i], ".txt")),
              sep="\t", quote=FALSE, row.names=FALSE)
}
saveRDS(mm_list, file.path(cfg$output$processed,
  paste0("module_membership_", wgcna_run_tag, ".rds")))

get_hub_lipids <- function(mm_df, top_n=10) {
  modules <- unique(mm_df$Module); modules <- modules[modules != "grey"]
  lapply(modules, function(mod) {
    sub <- mm_df[mm_df$Module == mod, ]
    if (!(mod %in% colnames(sub))) return(NULL)
    top_idx <- order(sub[[mod]], decreasing=TRUE)[1:min(top_n, nrow(sub))]
    data.frame(Module=mod, HubLipid=sub$Probe[top_idx],
               Rank=seq_along(top_idx), kME=sub[[mod]][top_idx],
               stringsAsFactors=FALSE)
  }) %>% bind_rows()
}

plot_hub_lipids <- function(hub_df, title_str, save_path) {
  hub_df$Module   <- fct_inorder(hub_df$Module)
  hub_df$HubLipid <- substr(hub_df$HubLipid, 1, 15)
  mod_labs    <- distinct(hub_df, Module)
  fill_colors <- setNames(levels(hub_df$Module), levels(hub_df$Module))
  text_colors <- fill_colors
  light_mods  <- c("lightcyan","green","cyan","lightyellow","yellow",
                   "grey60","pink","tan","lightgreen","salmon","greenyellow")
  text_colors[intersect(names(text_colors), light_mods)] <- "black"
  p <- ggplot(hub_df, aes(x=Rank, y=Module)) +
    geom_tile(data=mod_labs, aes(x=0.5, y=Module, fill=Module),
              width=0.4, height=0.8, inherit.aes=FALSE) +
    geom_text(aes(label=HubLipid, color=Module), size=3.5, fontface="bold") +
    scale_fill_manual(values=fill_colors) +
    scale_color_manual(values=text_colors) +
    scale_x_continuous(breaks=1:10, name="Rank", expand=expansion(add=c(0.5,0.5))) +
    labs(title=title_str) + theme_minimal() +
    theme(axis.text.y=element_text(face="bold"),
          panel.grid=element_blank(), legend.position="none")
  ggsave(save_path, plot=p, width=19.19, height=11.26, units="in")
}

title_map <- c(Baseline="Top 10 Hub Lipids (10-16 weeks)",
               Wk36_38="Top 10 Hub Lipids (36-38 weeks)",
               Postpartum="Top 10 Hub Lipids (3 Months Postpartum)")
for (tp in tp_labels) {
  hubs <- get_hub_lipids(mm_list[[tp]])
  write.csv(hubs, file.path(out, paste0("hub_lipids_top10_", tp, ".csv")), row.names=FALSE)
  plot_hub_lipids(hubs, paste0(title_map[[tp]], " [", wgcna_run_tag, "]"),
                  file.path(out, paste0("hub_lipids_", tp, ".pdf")))
}
message("Script 06 complete → ", out)