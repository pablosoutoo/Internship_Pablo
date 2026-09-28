set.seed(42)
#Module 1 - s2: annotation of the tumour-cell clusters (Harmony and CCA)
#Input : the integrated object saved at the end of s1_script.R
#Goal  : give every cluster a biological label, using
#        (1) published melanoma tumour-state signatures (scored per cell, averaged per cluster)
#        (2) canonical marker genes and data-driven markers (FindAllMarkers)
#The automatic label is a PROPOSAL: check the heatmap, dot plot and marker tables,
#and correct it in the "Manual review" section before using the labels downstream.

library(dplyr)
library(Seurat)
library(qs2)
library(ggplot2)

in_obj   <- "Module1 (Exploratory Analysis)/intermediate/s1_integrated_harmony_cca.qs2"
sig_file <- "Inputs/External/Tumor_Signatures.csv"   #PeeperLab/HeterotypicClustersR: extdata/tumor_data/Tumor_Signatures.csv
res_dir  <- "Module1 (Exploratory Analysis)/results/Annotation"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

data <- qs_read(in_obj)
DefaultAssay(data) <- "SCT"

#1. Signature scores per cell
##Tumor_Signatures.csv contains several published signature sets.
##"Pozniak" (Pozniak et al., Cell 2024, the TCF4 paper in the reading list) defines the tumour states
##used in the lab's own annotation (Melanocytic, Mitotic, Neural crest-like, Mesenchymal-like,
##IFN response, Antigen presentation, Stress hypoxia / p53...), so it is the main reference.
##"Tsoi" (differentiation stages along the melanocytic axis) is added as a second opinion.
sigs <- read.csv(sig_file)
sigs <- sigs[grepl("^(Pozniak|Tsoi) ", sigs$Signature), ]
sig_list <- split(sigs$Gene, sigs$Signature)
sig_list <- lapply(sig_list, function(g) intersect(unique(g), rownames(data)))
sapply(sig_list, length)                            #genes found in our data per signature
sig_list <- sig_list[sapply(sig_list, length) >= 10]
names(sig_list) <- make.names(names(sig_list))      #safe metadata column names

##AddModuleScore = mean expression of the signature minus the mean of random control genes
##with similar expression, so scores are comparable between cells.
data <- AddModuleScore(data, features = sig_list, name = "sig_", assay = "SCT", seed = 42)
score_cols <- paste0("sig_", seq_along(sig_list))
colnames(data@meta.data)[match(score_cols, colnames(data@meta.data))] <- names(sig_list)

pretty_name <- function(x) trimws(gsub("\\.+", " ", sub("^Pozniak\\.", "", x)))

#2. Automatic label per cluster
##Signatures have different sizes, so their raw scores are not comparable with each other.
##We therefore z-score each signature ACROSS clusters (which clusters are most enriched for it)
##and give each cluster the signature in which it stands out the most.
##Patient_specific_A/B are left out of the label pool: they describe patients of the Pozniak
##cohort, not ours. The margin (top z minus second z) says how clear-cut the label is.
label_pool <- grep("^Pozniak", names(sig_list), value = TRUE)
label_pool <- label_pool[!grepl("Patient_specific", label_pool)]

cluster_means <- function(cluster_col, cols) {
  data@meta.data %>%
    group_by(cluster = .data[[cluster_col]]) %>%
    summarise(across(all_of(cols), mean), n_cells = n(), .groups = "drop")
}

label_clusters <- function(cluster_col) {
  means <- cluster_means(cluster_col, label_pool)
  z <- scale(as.matrix(means[, label_pool]))
  ord <- apply(z, 1, function(v) names(sort(v, decreasing = TRUE)))
  data.frame(cluster = as.character(means$cluster),
             n_cells = means$n_cells,
             label   = pretty_name(ord[1, ]),
             second  = pretty_name(ord[2, ]),
             margin  = round(apply(z, 1, function(v) diff(sort(v, decreasing = TRUE)[2:1])), 2))
}

ann_harmony <- label_clusters("clusters_harmony")
ann_cca     <- label_clusters("clusters_cca")
ann_harmony
ann_cca

#3. Manual review
##Look at the heatmap, the dot plot and the marker tables below, then correct labels here.
##Example: manual_harmony <- c(`7` = "Mitotic", `11` = "Patient specific (P3)")
manual_harmony <- character(0)
manual_cca     <- character(0)
ann_harmony$label[match(names(manual_harmony), ann_harmony$cluster)] <- manual_harmony
ann_cca$label[match(names(manual_cca), ann_cca$cluster)] <- manual_cca

write.csv(ann_harmony, file.path(res_dir, "s2_cluster_labels_harmony.csv"), row.names = FALSE)
write.csv(ann_cca,     file.path(res_dir, "s2_cluster_labels_cca.csv"),     row.names = FALSE)

##Add the labels to every cell ("<cluster>: <state>" keeps clusters with the same state apart)
data$state_harmony <- ann_harmony$label[match(as.character(data$clusters_harmony), ann_harmony$cluster)]
data$state_cca     <- ann_cca$label[match(as.character(data$clusters_cca), ann_cca$cluster)]
data$cluster_label_harmony <- paste0(data$clusters_harmony, ": ", data$state_harmony)
data$cluster_label_cca     <- paste0(data$clusters_cca, ": ", data$state_cca)

#4. Evidence plots
##4a. Heatmap: mean signature score per cluster, z-scored across clusters (Pozniak + Tsoi)
heat <- function(cluster_col, title) {
  means <- cluster_means(cluster_col, names(sig_list))
  zmat  <- scale(as.matrix(means[, names(sig_list)]))
  long  <- data.frame(cluster   = rep(means$cluster, times = ncol(zmat)),
                      signature = rep(gsub("\\.+", " ", colnames(zmat)), each = nrow(zmat)),
                      z         = as.vector(zmat))
  ggplot(long, aes(cluster, signature, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "darkred",
                         name = "Mean score\n(z across clusters)") +
    labs(x = "Cluster", y = NULL, title = title) +
    theme_classic()
}
heat_h <- heat("clusters_harmony", "Harmony clusters - signature enrichment")
heat_c <- heat("clusters_cca", "CCA clusters - signature enrichment")
heat_h
heat_c
ggsave(file.path(res_dir, "s2_signature_heatmap_harmony.png"), heat_h, width = 10, height = 8, dpi = 300)
ggsave(file.path(res_dir, "s2_signature_heatmap_cca.png"),     heat_c, width = 10, height = 8, dpi = 300)

##4b. Dot plot of canonical genes for each state
canonical <- c("MKI67", "TOP2A",                          #cycling
               "MITF", "PMEL", "DCT", "MLANA", "TYR",     #melanocytic
               "SOX10", "NGFR", "AXL", "SOX9",            #neural crest-like / dedifferentiated
               "VIM", "SERPINE1", "FN1",                  #mesenchymal
               "VEGFA", "CA9", "NDUFA4L2",                #hypoxia
               "CXCL10", "STAT1", "B2M", "HLA-A",         #IFN response / MHC-I
               "CD74", "HLA-DRA",                         #MHC-II
               "HSPA6", "HSPA1A", "CDKN1A", "GDF15")      #stress (heat shock, p53)
canonical <- intersect(canonical, rownames(data))
dot_h <- DotPlot(data, features = canonical, group.by = "clusters_harmony", assay = "SCT") +
  RotatedAxis() + ggtitle("Harmony clusters - canonical genes")
dot_c <- DotPlot(data, features = canonical, group.by = "clusters_cca", assay = "SCT") +
  RotatedAxis() + ggtitle("CCA clusters - canonical genes")
dot_h
dot_c
ggsave(file.path(res_dir, "s2_canonical_dotplot_harmony.png"), dot_h, width = 12, height = 6, dpi = 300)
ggsave(file.path(res_dir, "s2_canonical_dotplot_cca.png"),     dot_c, width = 12, height = 6, dpi = 300)

##4c. Data-driven markers (top 20 per cluster). PrepSCTFindMarkers is needed because there is one SCT
##model per patient. max.cells.per.ident keeps the run time reasonable on Hooke.
data <- PrepSCTFindMarkers(data)
Idents(data) <- "clusters_harmony"
mk_h <- FindAllMarkers(data, assay = "SCT", only.pos = TRUE, min.pct = 0.25,
                       logfc.threshold = 0.5, max.cells.per.ident = 1000)
Idents(data) <- "clusters_cca"
mk_c <- FindAllMarkers(data, assay = "SCT", only.pos = TRUE, min.pct = 0.25,
                       logfc.threshold = 0.5, max.cells.per.ident = 1000)
write.csv(mk_h %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 20),
          file.path(res_dir, "s2_top20_markers_harmony.csv"), row.names = FALSE)
write.csv(mk_c %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 20),
          file.path(res_dir, "s2_top20_markers_cca.csv"), row.names = FALSE)

#5. Labelled UMAPs
umap_h <- DimPlot(data, reduction = "umap.harmony", group.by = "cluster_label_harmony",
                  label = TRUE, repel = TRUE, label.size = 3) + ggtitle("Harmony - proposed tumour states")
umap_c <- DimPlot(data, reduction = "umap.cca", group.by = "cluster_label_cca",
                  label = TRUE, repel = TRUE, label.size = 3) + ggtitle("CCA - proposed tumour states")
umap_h
umap_c
ggsave(file.path(res_dir, "s2_UMAP_states_harmony.png"), umap_h, width = 12, height = 8, dpi = 300)
ggsave(file.path(res_dir, "s2_UMAP_states_cca.png"),     umap_c, width = 12, height = 8, dpi = 300)

#6. Do the two methods give the same label to the same cell?
label_agreement <- mean(data$state_harmony == data$state_cca)
label_agreement
write.csv(as.data.frame.matrix(table(Harmony = data$state_harmony, CCA = data$state_cca)),
          file.path(res_dir, "s2_label_agreement_harmony_vs_cca.csv"))

#7. Save (large files: intermediate/ is git-ignored, sync to SurfDrive)
write.csv(data@meta.data[, c("patient_id", "sample_id", "clusters_harmony", "state_harmony",
                             "clusters_cca", "state_cca")],
          "Module1 (Exploratory Analysis)/intermediate/s2_cell_annotations.csv")
qs_save(data, "Module1 (Exploratory Analysis)/intermediate/s2_annotated.qs2")
