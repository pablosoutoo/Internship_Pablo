set.seed(42)
#Download the necessary libraries
library(dplyr)
library(Seurat)
library(SeuratData)
library(qs2)
library(ggplot2)
library(sctransform)

#Read Data
data<-qs_read("Inputs/External/stripped_s1_harmony_tumor_cells.qs2")
data

#Standard preprocessing workflow

##Calculate the mitochondrial QC metrics 
data[["percent.mt"]] <-PercentageFeatureSet(data, pattern = "^MT-")

##Visualize in plots
violin_plot<-VlnPlot(data, features =c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0, ncol=3)

ggsave(filename="Module1 (Exploratory Analysis)/results/s1_mitocond_pres.png", plot=violin_plot,width = 8, height = 6, dpi = 300)

##FeatureScatter <- Tipically used to visualize feature-feature relationship

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")

ggsave(filename="Module1 (Exploratory Analysis)/results/s1_feat_scatt_1.png", plot=plot1,width = 8, height = 6, dpi = 300)
ggsave(filename="Module1 (Exploratory Analysis)/results/s1_feat_scatt_2.png", plot=plot2,width = 8, height = 6, dpi = 300)

#We filter the data in order to discard cells with a low gene expression and a high percentage of mithocondrial genes.
#We'll keep the cells with at leats 500 genes detected and maximal 7500(nFeatureRNA), and below a 15% percentage of mithocondrial genes (already in the paper). 

# NOTE (2026-09-23): subset() is unreliable on this object. It trims meta.data
# correctly (37,407 QC-passing cells) but leaves the RNA assay's counts matrix at
# the original 45,747 cells, so NormalizeData() fails because there is a dimension mismatch.
# Ruled out: the ADT assay's un-joined per-sample layers (dropping ADT did not fix it),
# a stale object (re-running subset() fresh did not fix it), and duplicated or
# mismatched barcodes (zero duplicates, perfect 1:1 match between meta.data and RNA).
# Most likely a version-specific bug or corrupted internal cell bookkeeping in the
# Assay5 object (separate from colnames(), which is fine).
# Workaround: skip subset(). Filter the counts matrix with plain column indexing
# and rebuild a fresh Seurat object with CreateSeuratObject().
# ADT is not used downstream, so it is still dropped.

data[["ADT"]] <- NULL

#Confirm that we have only RNA left
Assays(data)
data[["RNA"]]

#For the record I'll write down the versions and assay class involved in the subset() bug
class(data[["RNA"]])
packageVersion("Seurat")
packageVersion("SeuratObject")

##QC-passing barcodes, taken from meta.data
qc_pass <- with(data@meta.data, nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 15)
keep_cells <- rownames(data@meta.data)[qc_pass]
length(keep_cells)

##Extract the raw counts (join per-sample layers first if the RNA assay is split)
if (length(Layers(data[["RNA"]], search = "counts")) > 1) {
  data[["RNA"]] <- JoinLayers(data[["RNA"]])
}
counts <- LayerData(data, assay = "RNA", layer = "counts")
stopifnot(!anyDuplicated(colnames(counts)), all(keep_cells %in% colnames(counts)))

##Subset by barcode with plain matrix indexing
counts_filtered <- counts[, keep_cells]

##Carry over metadata for the kept cells (nCount/nFeature are recomputed by CreateSeuratObject)
meta_filtered <- data@meta.data[keep_cells, setdiff(colnames(data@meta.data), c("nCount_RNA", "nFeature_RNA")), drop = FALSE]

##Rebuild a clean, RNA-only object (min.cells/min.features = 0 so no extra filtering happens)
data <- CreateSeuratObject(counts = counts_filtered, meta.data = meta_filtered,
                                    assay = "RNA", min.cells = 0, min.features = 0)
rm(counts, counts_filtered, meta_filtered)


##Sanity check: all three must agree (28,126 genes x 44,478 cells)
ncol(data)
nrow(data@meta.data)
dim(data[["RNA"]])
stopifnot(
  ncol(data) == length(keep_cells),
  nrow(data@meta.data) == length(keep_cells),
  ncol(data[["RNA"]]) == length(keep_cells),
  identical(colnames(data), rownames(data@meta.data))
)

## Mark doublets tumor–T (the object still has only one layer counts)
cnt   <- LayerData(data, assay = "RNA", layer = "counts")
t_id  <- c("CD3D", "CD3E", "CD2", "CD8A")
cyto  <- c("NKG7", "GZMB", "GZMA", "CCL5", "PRF1", "CST7")
data$T_doublet <- colSums(cnt[t_id, ] > 0) >= 2 & colSums(cnt[cyto, ] > 0) >= 2
table(data$T_doublet, data$sample_id)
rm(cnt)
data <- subset(data, subset = T_doublet == FALSE)   # works: object already built
data[["RNA"]] <- split(data[["RNA"]], f = data$patient_id)
stopifnot(length(unique(data$patient_id)) > 1)

#Normalization of the data with SCTransform()
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = 'percent.mt')
data<- SCTransform(data)

drop <- c(grep("^(MT-|RP[SL]\\d|TR[ABDG][VJC]|IG[HKL][VJC])", rownames(data), value = TRUE),
          "HTN1", "HTN3", "STATH", t_id, cyto,
          "PTPRC", "IL32", "CD52", "CORO1A", "CCL4", "CD7", "LAG3", "CD69")
VariableFeatures(data) <- setdiff(VariableFeatures(data), drop)

dim(data)

#Identification of highly variable features

#PCA
data <- RunPCA(data, features = VariableFeatures(object = data))
pca1<-VizDimLoadings(data,dim =1:2, reduction ="pca")
pca1
ggsave(filename="Module1 (Exploratory Analysis)/results/PCA/most_imp_genes.png", plot=pca1,width = 8, height = 6, dpi = 300)

pca2<-DimPlot(data, reduction = "pca")
pca2
ggsave(filename="Module1 (Exploratory Analysis)/results/PCA/PCA.png", plot=pca2,width = 8, height = 6, dpi = 300)

DimHeatmap(data, dims = 1, cells = 500, balanced = TRUE)

elbow_plot<-ElbowPlot(data)
elbow_plot

ggsave(filename="Module1 (Exploratory Analysis)/results/PCA/Elbow_plot.png", plot=elbow_plot,width = 8, height = 6, dpi = 300)

#Integration
#We integrate the same object with two different methods so we can compare them afterwards.
#Both start from the same SCT-normalised, per-patient layers and the same PCA, so any
#difference between the results comes from the integration method itself.

##Harmony (cluster-level correction of the PCA embedding)
data <- IntegrateLayers(object = data, method = HarmonyIntegration,
                        orig.reduction = "pca", new.reduction = "harmony",
                        normalization.method = "SCT", verbose = FALSE)

harmony_integration <-function(data,theta,lambda,max_iter,sigma){
   
   data <- IntegrateLayers(
      object = data, 
      method = HarmonyIntegration,
      orig.reduction = "pca",
      new.reduction = "harmony",
      normalization.method = "SCT",
      verbose = TRUE,
      
      #---Custom Harmony Parameters ---
      theta = theta,
      lambda = lambda,
      max.iter.harmony = max_iter,
      sigma = sigma
      
    )
   
    data <- FindNeighbors(data, reduction = "harmony", dims = 1:15,
                          graph.name = c("harmony_nn", "harmony_snn"))
    data <- FindClusters(data, graph.name = "harmony_snn", resolution = 0.5,
                         cluster.name = "clusters_harmony")
    
    data <- RunUMAP(data, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony")
    
    #A text to write down the parameters values
    params_text <- paste0("theta = ",theta, " | lambda = ", lambda, " | max.iter.harmony =", max_iter, " | sigma = ", sigma)
     
    #A filename with the values
    file_name <- paste0("Module1 (Exploratory Analysis)/results/UMAP/UMAP_theta", theta,
                        "_lambda", lambda, "_max.iter.harmony" ,max_iter, "_sigma",sigma  , ".png")
    
    
    umap_plot <- DimPlot(data, reduction = "umap.harmony",
                        group.by = "clusters_harmony", label = TRUE) +
      ggtitle(params_text)
    
    ggsave(filename = file_name, plot = umap_plot, width = 8, height = 6, dpi = 300)
    
    return(umap_plot)
}

test_plot <- harmony_integration(data, theta = 2, lambda = 1, max_iter = 10, sigma = 0.1)
test_plot


##CCA (Seurat anchors: cell-level mutual nearest neighbours in a shared CCA space)
##This is the method used in the original tumour analysis (Tumor_Analysis_Code.R).
##It is slower and uses more memory than Harmony, so it may take a while on Hooke.
data <- IntegrateLayers(object = data, method = CCAIntegration,
                        orig.reduction = "pca", new.reduction = "integrated.cca",
                        normalization.method = "SCT", verbose = FALSE)



##Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data),10)
top10

##Plot variable features with and without labels
hvf <- SCTResults(data[["SCT"]], slot = "feature.attributes")
if (is.data.frame(hvf)) hvf <- list(all = hvf)
hvf_df <- do.call(rbind, lapply(names(hvf), function(m) data.frame(
  model = m, gene = rownames(hvf[[m]]),
  gmean = hvf[[m]]$gmean, residual_variance = hvf[[m]]$residual_variance)))
hvf_df$variable <- hvf_df$gene %in% VariableFeatures(data)
plot2a <- ggplot(hvf_df, aes(gmean, residual_variance, colour = variable)) +
  geom_point(size = 0.3) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "red"),
                      labels = c("no", "yes"), name = "Variable") +
  ggrepel::geom_text_repel(data = subset(hvf_df, gene %in% top10), aes(label = gene),
                           colour = "black", size = 2.5, max.overlaps = Inf) +
  facet_wrap(~ model) +
  labs(x = "Geometric mean of expression", y = "Residual variance") +
  theme_classic()
plot2a

ggsave(filename="Module1 (Exploratory Analysis)/results/s1_variable_features.png", plot=plot2a,width = 8, height = 6, dpi = 300)

#Clusterization
#Same dims and resolution for both methods so the comparison is fair.
#Each method gets its own graph and cluster column, so nothing is overwritten.
#(Do not use "harmony_clusters" as a name: that column holds the lab's original clusters.)

##Harmony
data <- FindNeighbors(data, reduction = "harmony", dims = 1:15,
                      graph.name = c("harmony_nn", "harmony_snn"))
data <- FindClusters(data, graph.name = "harmony_snn", resolution = 0.5,
                     cluster.name = "clusters_harmony")

##CCA
data <- FindNeighbors(data, reduction = "integrated.cca", dims = 1:15,
                      graph.name = c("cca_nn", "cca_snn"))
data <- FindClusters(data, graph.name = "cca_snn", resolution = 0.5,
                     cluster.name = "clusters_cca")

table(data$clusters_harmony)
table(data$clusters_cca)

#UMAP/t-SNE
data <- RunUMAP(data, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony")
data <- RunUMAP(data, reduction = "integrated.cca", dims = 1:15, reduction.name = "umap.cca")

umap_plot<-DimPlot(data, reduction = "umap.harmony", group.by = "patient_id", shuffle = TRUE)
umap_plot
ggsave(filename="Module1 (Exploratory Analysis)/results/UMAP/UMAP.png", plot=umap_plot,width = 8, height = 6, dpi = 300)

#Comparison of the two integration methods
comp_dir <- "Module1 (Exploratory Analysis)/results/Integration_comparison"
dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

##1. UMAPs side by side: top row coloured by patient (batch mixing), bottom row by cluster
umap_compare <- patchwork::wrap_plots(
  DimPlot(data, reduction = "umap.harmony", group.by = "patient_id", shuffle = TRUE) + ggtitle("Harmony - patient"),
  DimPlot(data, reduction = "umap.cca", group.by = "patient_id", shuffle = TRUE) + ggtitle("CCA - patient"),
  DimPlot(data, reduction = "umap.harmony", group.by = "clusters_harmony", label = TRUE) + NoLegend() + ggtitle("Harmony - clusters"),
  DimPlot(data, reduction = "umap.cca", group.by = "clusters_cca", label = TRUE) + NoLegend() + ggtitle("CCA - clusters"),
  ncol = 2)
umap_compare
ggsave(filename = file.path(comp_dir, "s1_UMAP_harmony_vs_cca.png"), plot = umap_compare, width = 14, height = 12, dpi = 300)

##2. Agreement between clusterings: Adjusted Rand Index (1 = identical, ~0 = random)
##Written in base R so no extra package is needed on Hooke.
ari <- function(x, y) {
  tab      <- table(x, y)
  sum_ij   <- sum(choose(tab, 2))
  sum_a    <- sum(choose(rowSums(tab), 2))
  sum_b    <- sum(choose(colSums(tab), 2))
  expected <- sum_a * sum_b / choose(sum(tab), 2)
  (sum_ij - expected) / ((sum_a + sum_b) / 2 - expected)
}
ari_table <- data.frame(
  comparison = c("Harmony vs CCA",
                 "Harmony vs original lab clusters",
                 "CCA vs original lab clusters"),
  ARI = c(ari(data$clusters_harmony, data$clusters_cca),
          ari(data$clusters_harmony, data$harmony_clusters),
          ari(data$clusters_cca,     data$harmony_clusters)))
ari_table
write.csv(ari_table, file.path(comp_dir, "s1_ARI_integration_methods.csv"), row.names = FALSE)

##3. Which Harmony cluster corresponds to which CCA cluster (row-normalised)
cross <- as.data.frame(prop.table(table(harmony = data$clusters_harmony,
                                        cca = data$clusters_cca), margin = 1))
cross_plot <- ggplot(cross, aes(cca, harmony, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred", name = "Fraction of\nHarmony cluster") +
  labs(x = "CCA cluster", y = "Harmony cluster") +
  theme_classic()
cross_plot
ggsave(filename = file.path(comp_dir, "s1_cluster_crosstab_harmony_vs_cca.png"), plot = cross_plot, width = 8, height = 7, dpi = 300)

##4. Patient mixing per cluster: normalised Shannon entropy of patient_id
##(1 = cells evenly spread over all patients, 0 = cluster made of a single patient).
##Clusters close to 0 are patient-specific; compare how many each method keeps.
patient_mixing <- function(clusters, patients, method) {
  tab <- table(clusters, patients)
  p   <- prop.table(tab, margin = 1)
  ent <- apply(p, 1, function(x) { x <- x[x > 0]; -sum(x * log(x)) }) / log(ncol(tab))
  data.frame(method = method, cluster = rownames(tab), n_cells = rowSums(tab),
             mixing = ent, top_patient = colnames(tab)[apply(tab, 1, which.max)])
}
mixing <- rbind(patient_mixing(data$clusters_harmony, data$patient_id, "Harmony"),
                patient_mixing(data$clusters_cca,     data$patient_id, "CCA"))
mixing
write.csv(mixing, file.path(comp_dir, "s1_patient_mixing_per_cluster.csv"), row.names = FALSE)

##5. Patient composition of each cluster
comp_plot <- function(cluster_col, title) {
  ggplot(data@meta.data, aes(x = .data[[cluster_col]], fill = patient_id)) +
    geom_bar(position = "fill") +
    labs(x = "Cluster", y = "Fraction of cells", title = title) +
    theme_classic()
}
composition <- patchwork::wrap_plots(comp_plot("clusters_harmony", "Harmony"),
                                     comp_plot("clusters_cca", "CCA"),
                                     ncol = 1, guides = "collect")
composition
ggsave(filename = file.path(comp_dir, "s1_patient_composition_per_cluster.png"), plot = composition, width = 10, height = 8, dpi = 300)

#Save the integrated object so s2_script.R (annotation) can start from here without re-running s1.
#It is a large file: intermediate/ is git-ignored, sync it to SurfDrive with rclone instead.
dir.create("Module1 (Exploratory Analysis)/intermediate", showWarnings = FALSE)
qs_save(data, "Module1 (Exploratory Analysis)/intermediate/s1_integrated_harmony_cca.qs2")
