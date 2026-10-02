set.seed(42)
#Module 1 - s2: annotation of the tumour-cell clusters (Harmony and CCA)
#Input : the integrated object saved at the end of s1_script.R
#Goal  : give every cluster a biological label, using
#        (1) published melanoma tumour-state signatures (scored per cell, averaged per cluster)
#        (2) canonical marker genes and data-driven markers (AUC-ranked and conserved across patients)
#        and consolidate the clusters into 7-8 meta-states for the figures
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
##The CCA labels are the ones first set by hand in s1_script.R (canonical dot plot + top markers).
manual_harmony <- character(0)
manual_cca     <- c(`0`  = "Neural crest-like",             #SOX10 high; NCMAP, SCN7A, SCRG1, ANGPTL7
                    `1`  = "Melanocytic",                   #PMEL, MLANA, MITF, DCT, TYR
                    `2`  = "Stress (ATF4 / amino acid)",    #ASNS, TRIB3, GDF15, CDKN1A (weak markers)
                    `3`  = "IFN response",                  #GBP1/4, IFIT2, IFI44L, STAT1, B2M, HLA-A
                    `4`  = "Mesenchymal-like (invasive)",   #MMP1, MMP3, IL11, SERPINB2, SERPINE1, INHBA
                    `5`  = "Stress (hypoxia)",              #NDUFA4L2, VEGFA, CA9, MT3
                    `6`  = "Melanocytic (pigmentation)",    #TYR, MITF, DCT high; low MHC-I
                    `7`  = "Mesenchymal-like (TGFb/YAP)",   #FN1, TAGLN, CCN1, CCN2, DKK1
                    `8`  = "Mitotic (G1/S)",                #E2F2, RRM2, MCM10, CDC45, CLSPN
                    `9`  = "Antigen presentation (MHC-II)", #CD74, HLA-DRA; B-cell genes in ~5% of cells
                    `10` = "Inflammatory (NF-kB)",          #CXCL10/11, CCL2, CXCL2, SELE, HSPA6
                    `11` = "Mitotic (G2/M)")                #PLK1, CDC20, KIF20A, MKI67, TOP2A
manual_cca <- manual_cca[names(manual_cca) %in% ann_cca$cluster]
ann_harmony$auto_label <- ann_harmony$label   #keep the signature-based label to compare
ann_cca$auto_label     <- ann_cca$label
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

##4c. Data-driven markers (top 20 per cluster), ranked by AUC instead of avg_log2FC.
##Ranking by fold change alone pushed up genes seen in very few cells (e.g. NCMAP in 1.8% of
##cluster 0 in s1_top10_markers_cca.csv). AUC rewards genes that actually separate the cluster.
##PrepSCTFindMarkers is needed because there is one SCT model per patient (without it,
##FindAllMarkers fails silently for every cluster and returns an EMPTY table).
##presto (immunogenomics/presto) computes the Wilcoxon test + AUC for all clusters in seconds.
##Note: presto's logFC is a difference of means of log-normalised data (natural log), not log2.
data <- PrepSCTFindMarkers(data)
expr <- GetAssayData(data, assay = "SCT", layer = "data")

auc_markers <- function(cluster_col) {
  presto::wilcoxauc(expr, as.character(data@meta.data[[cluster_col]])) %>%
    filter(padj < 0.05, logFC > 0.25, pct_in >= 25, auc >= 0.6) %>%
    mutate(pct_diff = pct_in - pct_out) %>%
    rename(cluster = group, gene = feature) %>%
    group_by(cluster) %>%
    arrange(desc(auc), .by_group = TRUE)
}
mk_h <- auc_markers("clusters_harmony")
mk_c <- auc_markers("clusters_cca")
rm(expr)
stopifnot(nrow(mk_h) > 0, nrow(mk_c) > 0)
mk_c %>% count(cluster)   #clusters with few or no markers passing the filters are suspicious
write.csv(mk_h %>% slice_head(n = 20), file.path(res_dir, "s2_top20_markers_auc_harmony.csv"), row.names = FALSE)
write.csv(mk_c %>% slice_head(n = 20), file.path(res_dir, "s2_top20_markers_auc_cca.csv"),     row.names = FALSE)

##4d. Markers conserved across patients (CCA clusters). A gene is kept only if it is up in the
##cluster in EVERY patient that has cells there (min_log2FC > 0); n_patients = 1 means the
##cluster is patient-specific. Needs the "metap" package. Patients with < 3 cells in a
##cluster are skipped with a warning. If you get "multiple models with unequal library sizes",
##add recorrect_umi = FALSE.
Idents(data) <- "clusters_cca"
conserved <- lapply(levels(data$clusters_cca), function(cl) {
  m <- tryCatch(FindConservedMarkers(data, ident.1 = cl, grouping.var = "patient_id",
                                     assay = "SCT", only.pos = TRUE, min.pct = 0.25,
                                     logfc.threshold = 0.25, max.cells.per.ident = 500,
                                     verbose = FALSE),
                error = function(e) { message("cluster ", cl, ": ", conditionMessage(e)); NULL })
  if (is.null(m) || nrow(m) == 0) return(NULL)
  fc_cols <- grep("_avg_log2FC$", colnames(m), value = TRUE)
  m$min_log2FC <- apply(m[, fc_cols, drop = FALSE], 1, min)
  m$n_patients <- length(fc_cols)
  data.frame(cluster = cl, gene = rownames(m),
             m[, intersect(c("min_log2FC", "n_patients", "max_pval", "minimump_p_val"), colnames(m))],
             row.names = NULL) %>%
    filter(min_log2FC > 0) %>%
    arrange(desc(min_log2FC)) %>%
    head(20)
})
conserved <- do.call(rbind, conserved)
conserved
write.csv(conserved, file.path(res_dir, "s2_top20_conserved_markers_cca.csv"), row.names = FALSE)

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

#7. Meta-states: consolidate the 12 CCA clusters into 7-8 tumour states (Pozniak et al. 2024)
##PROVISIONAL. Decide with: the signature heatmap and auto_label (sections 2 and 4a), the AUC and
##conserved markers (4c, 4d) and the QC by cluster from s1 (results/Diagnostics/s1_qc_by_cluster.csv).
##   cluster 2  -> "Low quality" if nFeature is low / percent.mt high
##   cluster 9  -> "Excluded (doublets)" if B_frac / myeloid_frac is high
##   cluster 10 -> "Excluded (dissociation)" if dissoc_score1 is high
##Clusters labelled "Excluded..." or "Low quality" are left out of the meeting figures.
meta_map <- c(`1`  = "Melanocytic",          `6`  = "Melanocytic",
              `0`  = "Neural crest-like",
              `4`  = "Mesenchymal-like",     `7`  = "Mesenchymal-like",
              `8`  = "Mitotic",              `11` = "Mitotic",
              `3`  = "IFN response",
              `9`  = "Antigen presentation",
              `5`  = "Stress (hypoxia)",
              `2`  = "Stress (ATF4)",
              `10` = "Inflammatory (NF-kB)")
stopifnot(all(levels(data$clusters_cca) %in% names(meta_map)))
data$meta_state <- factor(unname(meta_map[as.character(data$clusters_cca)]),
                          levels = unique(meta_map))
table(data$clusters_cca, data$meta_state)

##Does the manual meta-state agree with the automatic Pozniak label of each cluster?
meta_check <- merge(ann_cca[, c("cluster", "n_cells", "auto_label", "second", "margin", "label")],
                    data.frame(cluster = names(meta_map), meta_state = unname(meta_map)), by = "cluster")
meta_check[order(as.numeric(meta_check$cluster)), ]
write.csv(meta_check, file.path(res_dir, "s2_meta_states_cca.csv"), row.names = FALSE)


#8. Figures for the meeting (UMAP, reduced dot plot, composition per patient)
##One fixed palette, so the same state has the same colour in every figure.
fig_dir <- "Module1 (Exploratory Analysis)/results/Figures_meeting"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

##The excluded cells are filtered inside each plot (no copy of the multi-GB object).
fig_states <- grep("^(Excluded|Low quality)", levels(data$meta_state), value = TRUE, invert = TRUE)
fig_cells  <- colnames(data)[data$meta_state %in% fig_states]
pal <- setNames(scales::hue_pal()(length(fig_states)), fig_states)

##A. Final UMAP coloured by meta-state
figA <- DimPlot(data, reduction = "umap.cca", group.by = "meta_state", cells = fig_cells, cols = pal,
                label = TRUE, repel = TRUE, label.size = 4) +
  NoLegend() + ggtitle(NULL)

##B. Reduced dot plot: 2-3 genes per meta-state, in the same order as the states
dot_genes <- c("PMEL", "MLANA", "MITF",         #melanocytic
               "SOX10", "NGFR",                 #neural crest-like
               "SERPINE1", "FN1",               #mesenchymal-like
               "MKI67", "TOP2A",                #mitotic
               "CXCL10", "STAT1",               #IFN response
               "CD74", "HLA-DRA",               #antigen presentation
               "VEGFA", "CA9",                  #hypoxia
               "ASNS", "TRIB3",                 #ATF4 stress
               "CXCL2", "HSPA6")                #NF-kB / inflammatory
Idents(data) <- "meta_state"
figB <- DotPlot(data, features = intersect(dot_genes, rownames(data)),
                idents = fig_states, assay = "SCT") +
  RotatedAxis() + labs(x = NULL, y = NULL)

##C. Composition per patient
comp <- as.data.frame(table(patient = data$patient_id[data$meta_state %in% fig_states],
                            state   = droplevels(data$meta_state[data$meta_state %in% fig_states])))
figC <- ggplot(comp, aes(patient, Freq, fill = state)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = "Fraction of tumour cells", fill = NULL) +
  theme_classic() + RotatedAxis()

figA
figB
figC
figs <- list(A_UMAP_meta_states    = list(figA, 7, 6),
             B_dotplot_meta_states = list(figB, 10, 5),
             C_composition_patient = list(figC, 8, 5))
for (f in names(figs)) {
  ggsave(file.path(fig_dir, paste0(f, ".pdf")), figs[[f]][[1]], width = figs[[f]][[2]], height = figs[[f]][[3]])
  ggsave(file.path(fig_dir, paste0(f, ".png")), figs[[f]][[1]], width = figs[[f]][[2]], height = figs[[f]][[3]], dpi = 300)
}
rm(figs)
gc()


#9. Save (large file: s2_annotated.qs2 is in .gitignore, sync it to SurfDrive)
write.csv(data@meta.data[, c("patient_id", "sample_id", "clusters_harmony", "state_harmony",
                             "clusters_cca", "state_cca", "meta_state")],
          "Module1 (Exploratory Analysis)/intermediate/s2_cell_annotations.csv")
qs_save(data, "Module1 (Exploratory Analysis)/intermediate/s2_annotated.qs2")
