##FGSEA ANALYSIS

#Libraries
library(fgsea)
library(data.table)
library(ggplot2)
library(msigdbr)
library(qs2)
library(dplyr)
library(Seurat)
library(MAST)

#Paths and settings
in_obj           <- "Module1 (Exploratory Analysis)/intermediate/data_clusterized.qs2"
fgsea_dir        <- "Module1 (Exploratory Analysis)/results/FGSEA"
intermediate_dir <- "Module1 (Exploratory Analysis)/intermediate"
dir.create(fgsea_dir, showWarnings = FALSE, recursive = TRUE)

target_cluster   <- "2"
de_test          <- "MAST"      #as in the lab's script
run_all_clusters <- FALSE      #TRUE = repeat the analysis for every cluster and draw an NES heatmap to compare them
padj_cutoff      <- 0.05


#Data
##Download the whole C2:CP collection for Human (Homo Sapiens)
##subcategory = "CP" only returns the generic CP sets, not Reactome, KEGG, WikiPathways, PID, BioCarta...
##so we take all of C2 and keep every subcollection starting with "CP".
##msigdbr >= 10 renamed category -> collection (and gs_subcat -> gs_subcollection)
if (packageVersion("msigdbr") >= "10.0.0") {
  c2_human <- msigdbr(species = "Homo sapiens", collection = "C2")
  c2_subcol <- c2_human$gs_subcollection
} else {
  c2_human <- msigdbr(species = "Homo sapiens", category = "C2")
  c2_subcol <- c2_human$gs_subcat
}
c2_cp_human <- c2_human[grepl("^CP", c2_subcol), ]
table(c2_subcol[grepl("^CP", c2_subcol)])   #check which subcollections were kept

data <- qs_read(in_obj)
Idents(data) <- "clusters_cca"   #PrepSCTFindMarkers() was already run in s1_script.R

#Named list for fgsea
pathways <- lapply(split(c2_cp_human$gene_symbol, c2_cp_human$gs_name), unique)


#Functions
#1 DE of one cluster vs all the other clusters
de_vs_rest <- function(cluster_id) {
  de <- FindMarkers(data, ident.1 = cluster_id, ident.2 = NULL, assay = "SCT",
                    test.use = de_test, latent.vars = "patient_id",
                    logfc.threshold = 0, min.pct = 0.05, only.pos = FALSE)
  de$gene <- rownames(de)
  de
}

#2. Ranked list
##p = 0 gives -log10(p) = Inf, and those (the most significant genes) would be dropped by is.finite(),
##so they get the smallest non-zero p-value instead
make_ranking <- function(de) {
  p <- de$p_val
  p[p == 0] <- min(p[p > 0])
  stats <- sign(de$avg_log2FC) * -log10(p)
  names(stats) <- de$gene
  stats <- stats[is.finite(stats)]
  sort(stats, decreasing = TRUE)
}

#3.FGSEA
run_fgsea <- function(stats, cluster_id) {
  set.seed(42)
  res <- fgsea(pathways = pathways, stats = stats, minSize = 15, maxSize = 500)
  res <- as.data.frame(res)
  res <- res[!is.na(res$padj), ]
  res$cluster <- cluster_id
  ord <- order(res$padj, -abs(res$NES))
  res <- res[ord, ]
  rownames(res) <- NULL
  res
}

#4. Main pathways: collapsePathways() removes the significant pathways that are only significant
##because they overlap with a stronger one (e.g. many Reactome sets sharing the same genes)
add_main_pathways <- function(res, stats) {
  res$main_pathway <- FALSE
  sig <- as.data.table(res[res$padj < padj_cutoff, ])
  if (nrow(sig) == 0) return(res)
  collapsed <- collapsePathways(sig[order(pval)], pathways, stats)
  res$main_pathway <- res$pathway %in% collapsed$mainPathways
  res
}

#5. Barplot with the top up- and down-regulated pathway
plot_top_nes <- function(res, cluster_id, n = 10, subtitle = paste0("* padj < ", padj_cutoff)) {
  df  <- as.data.frame(res)[, c("pathway", "NES", "padj")]
  pool <- if (any(df$padj < padj_cutoff)) df[df$padj < padj_cutoff, ] else df
  top <- bind_rows(pool %>% filter(NES > 0) %>% slice_max(NES, n = n, with_ties = FALSE),
                   pool %>% filter(NES < 0) %>% slice_min(NES, n = n, with_ties = FALSE))
  top$label <- gsub("_", " ", top$pathway)
  ggplot(top, aes(NES, reorder(label, NES), fill = NES)) +
    geom_col() +
    geom_point(data = top[top$padj < padj_cutoff, ], shape = 8, size = 1.5, show.legend = FALSE) +
    scale_fill_gradient2(low = "#7cd3f7", mid = "white", high = "#825d9e") +
    scale_y_discrete(labels = scales::label_wrap(45)) +
    labs(x = "NES", y = NULL, title = paste("Cluster", cluster_id, "vs rest - MSigDB C2:CP"), subtitle = subtitle) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 6), panel.grid = element_blank())
}


#ANALYSIS OF CLUSTER 2
##1. DE and ranking
de_target    <- de_vs_rest(target_cluster)
stats_target <- make_ranking(de_target)
head(stats_target, 20)   #most up-regulated in cluster 2
tail(stats_target, 20)   #most down-regulated in cluster 2

##2 FGSEA
fgsea_target <- run_fgsea(stats_target, target_cluster)
fgsea_target <- add_main_pathways(fgsea_target, stats_target)

##Save the table (fwrite writes the leadingEdge list column as genes separated by "|")
fwrite(fgsea_target, file.path(fgsea_dir, paste0("fgsea_cluster", target_cluster, ".csv")))

##3 Leading edge
top_up <- head(fgsea_target[fgsea_target$NES > 0 & fgsea_target$padj < padj_cutoff, ], 10)
data.frame(pathway = top_up$pathway,
           NES     = round(top_up$NES, 2),
           leading_edge = sapply(top_up$leadingEdge, function(g) paste(head(g, 15), collapse = ", ")))

#Plots
##Top up /down pathways
p_top <- plot_top_nes(fgsea_target, target_cluster)
p_top
ggsave(file.path(fgsea_dir, paste0("top_nes_cluster", target_cluster, ".png")), p_top, width = 8, height = 7, dpi = 300)


#Enrichment plots of the top 3 up-regulated main pathways
top3 <- head(fgsea_target$pathway[fgsea_target$NES > 0 & fgsea_target$main_pathway], 3)
for (pw in top3) {
  res_pw <- fgsea_target[fgsea_target$pathway == pw, ]
  p_enr <- plotEnrichment(pathways[[pw]], stats_target) +
    labs(title = gsub("_", " ", pw),
         subtitle = paste0("Cluster ", target_cluster, " | NES = ", round(res_pw$NES, 2),
                           " | padj = ", signif(res_pw$padj, 2)))
  print(p_enr)
  ggsave(file.path(fgsea_dir, paste0("enrichment_cluster", target_cluster, "_", pw, ".png")), p_enr, width = 7, height = 4, dpi = 300)
}


#ALL CLUSTERS (optional)
if (run_all_clusters) {
  clusters <- levels(Idents(data))
  fgsea_all <- bind_rows(lapply(clusters, function(cl) {
    message("Cluster ", cl)
    run_fgsea(make_ranking(de_vs_rest(cl)), cl)
  }))
  fwrite(fgsea_all, file.path(fgsea_dir, "fgsea_all_clusters.csv"))

  ##Heatmap with the top 5 significant pathways (by |NES|) of each cluster
  top_pw <- fgsea_all %>%
    filter(padj < padj_cutoff) %>%
    group_by(cluster) %>%
    slice_max(abs(NES), n = 5, with_ties = FALSE) %>%
    pull(pathway) %>%
    unique()
  hm <- fgsea_all %>% filter(pathway %in% top_pw)
  hm$cluster <- factor(hm$cluster, levels = clusters)
  hm$label   <- gsub("_", " ", hm$pathway)

  p_hm <- ggplot(hm, aes(cluster, label, fill = NES)) +
    geom_tile(colour = "grey80") +
    geom_text(data = hm[hm$padj < padj_cutoff, ], label = "*", size = 3) +
    scale_fill_gradient2(low = "#7cd3f7", mid = "white", high = "#825d9e") +
    scale_y_discrete(labels = scales::label_wrap(50)) +
    labs(x = "Cluster", y = NULL, title = "NES per cluster vs rest - MSigDB C2:CP",
         subtitle = paste0("* padj < ", padj_cutoff)) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 6), panel.grid = element_blank())
  print(p_hm)
  ggsave(file.path(fgsea_dir, "heatmap_nes_all_clusters.png"), p_hm, width = 9, height = 10, dpi = 300)
}


#Save session information
Info_script <-sessionInfo()

#Save environment
save.image("Module1 (Exploratory Analysis)/intermediate/environment_fgsea.rdata")
