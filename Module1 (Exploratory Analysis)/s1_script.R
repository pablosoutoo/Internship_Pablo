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

diag_dir <- "Module1 (Exploratory Analysis)/results/Diagnostics"
dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)

##Cell counts per patient at every filtering step, written to a table instead of hard-coded
##in comments (the old comments disagreed: 37,407 vs 44,478 QC-passing cells).
count_cells <- function(obj, step) {
  n <- table(obj$patient_id)
  data.frame(step = step, patient_id = c(names(n), "ALL"), cells = c(as.vector(n), ncol(obj)))
}
cell_log <- count_cells(data, "1_input")

##How many patients and samples are really in the object (compare with the paper)
pat_samp <- table(patient = data$patient_id, sample = data$sample_id)
pat_samp
rowSums(pat_samp > 0)   #samples per patient: if > 1, consider splitting the layers by sample_id instead of patient_id
write.csv(as.data.frame.matrix(pat_samp), file.path(diag_dir, "s1_patient_x_sample_input.csv"))

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
# correctly but leaves the RNA assay's counts matrix at
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


##Sanity check: all three must agree (the real numbers are stored in cell_log)
ncol(data)
nrow(data@meta.data)
dim(data[["RNA"]])
stopifnot(
  ncol(data) == length(keep_cells),
  nrow(data@meta.data) == length(keep_cells),
  ncol(data[["RNA"]]) == length(keep_cells),
  identical(colnames(data), rownames(data@meta.data))
)
cell_log <- rbind(cell_log, count_cells(data, "2_QC_filter"))

## Mark doublets tumor–T (the object still has only one layer counts)
cnt   <- LayerData(data, assay = "RNA", layer = "counts")
t_id  <- c("CD3D", "CD3E", "CD2", "CD8A")
cyto  <- c("NKG7", "GZMB", "GZMA", "CCL5", "PRF1", "CST7")
data$T_doublet <- colSums(cnt[t_id, ] > 0) >= 2 & colSums(cnt[cyto, ] > 0) >= 2
table(data$T_doublet, data$sample_id)

##Same rule for B-cell and myeloid doublets. These are only FLAGGED, not removed:
##first check how they spread over the clusters (QC by cluster, after the integration).
##CD74/HLA-DR are not used on purpose: melanoma cells can express MHC-II themselves.
b_id  <- intersect(c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1"), rownames(cnt))
my_id <- intersect(c("LYZ", "CD14", "C1QA", "C1QB", "AIF1", "TYROBP"), rownames(cnt))
data$B_doublet       <- colSums(cnt[b_id, ]  > 0) >= 2
data$myeloid_doublet <- colSums(cnt[my_id, ] > 0) >= 2
table(B = data$B_doublet, myeloid = data$myeloid_doublet)

##scDblFinder simulates doublets within each sample. Also only flagged here.
##Limitation: the object holds tumour cells only, so it finds tumour-tumour doublets
##(between states or clones); tumour-immune doublets are better caught by the rules above.
set.seed(42)
sce <- SingleCellExperiment::SingleCellExperiment(list(counts = cnt))
sce <- scDblFinder::scDblFinder(sce, samples = data$sample_id, BPPARAM = BiocParallel::SerialParam())
data$dbl_score <- sce$scDblFinder.score
data$dbl_class <- as.character(sce$scDblFinder.class)
table(data$dbl_class, data$sample_id)
rm(cnt, sce)
gc()

data <- subset(data, subset = T_doublet == FALSE)   # works: object already built
cell_log <- rbind(cell_log, count_cells(data, "3_T_doublets_removed"))
cell_counts <- xtabs(cells ~ patient_id + step, data = cell_log)
cell_counts
write.csv(as.data.frame.matrix(cell_counts), file.path(diag_dir, "s1_cells_per_step.csv"))
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

##Scores to tell real states from technical artefacts (used in the QC by cluster below)
##Dissociation stress: genes induced by tissue digestion (van den Brink et al. 2017).
##Ambient RNA: highly abundant transcripts of OTHER cell types (plasma cells, red cells,
##salivary gland). Low levels in many tumour cells = ambient; they are the reason why
##IG and HTN/STATH genes had to be removed from the variable features above.
dissoc  <- c("FOS", "FOSB", "JUN", "JUNB", "ATF3", "EGR1", "IER2", "IER3", "DUSP1", "ZFP36",
             "NR4A1", "HSPA1A", "HSPA1B", "HSPA6", "DNAJB1")
ambient <- c("IGKC", "IGHG1", "IGLC2", "IGHA1", "JCHAIN", "HBB", "HBA1", "HTN1", "HTN3", "STATH")
data <- AddModuleScore(data, features = list(intersect(dissoc, rownames(data))),  name = "dissoc_score",  seed = 42)
data <- AddModuleScore(data, features = list(intersect(ambient, rownames(data))), name = "ambient_score", seed = 42)

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

#The grid functions below return only a light per-cell table (UMAP coordinates + clusters),
#never the Seurat object or ggplot objects, so memory does not grow with every run.
#All figures are rebuilt afterwards from that table.


#Output folders
harmony_dir      <- "Module1 (Exploratory Analysis)/results/UMAP/harmony_grid"
cca_dir          <- "Module1 (Exploratory Analysis)/results/CCA/cca_grid"
intermediate_dir <- "Module1 (Exploratory Analysis)/intermediate"
for (d in c(harmony_dir, cca_dir, intermediate_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)



##Rebuilds the "clusters + patients" figure from one run of the grid table.
##Mimics DimPlot: same point size, cluster labels at the median of each cluster,
##and patients drawn in shuffled order so no patient hides the others.
plot_umap_run <- function(df, axis_prefix) {
  pt_size <- min(1583 / nrow(df), 1)
  x_lab <- paste0(axis_prefix, "1")
  y_lab <- paste0(axis_prefix, "2")
  params_text <- df$params[1]
  
  centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ cluster, data = df, FUN = median)
  
  plot_clusters <- ggplot(df, aes(UMAP_1, UMAP_2, colour = cluster)) +
    geom_point(size = pt_size, shape = 16) +
    geom_text(data = centers, aes(label = cluster), colour = "black", size = 4) +
    labs(x = x_lab, y = y_lab) +
    theme_classic() +
    NoLegend() 
  
  set.seed(42)
  df_shuffled <- df[sample(nrow(df)), ]
  plot_patients <- ggplot(df_shuffled, aes(UMAP_1, UMAP_2, colour = patient_id)) +
    geom_point(size = pt_size, shape = 16) +
    labs(x = x_lab, y = y_lab, colour = "patient_id") +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    theme_classic() 
  
  patchwork::wrap_plots(plot_clusters, plot_patients) +
    patchwork::plot_annotation(title = df$params[1],
                               theme = theme(plot.title = element_text(size = 11)))
}

##Overview of a whole grid in ONE ggplot, faceted by the parameters (built from the table,
##so it cannot pick up a stale grid or plot_list). colour_by = "cluster" or "patient_id".
plot_grid_overview <- function(results, facet_formula, colour_by, title) {
  facet_vars <- all.vars(facet_formula)
  if (colour_by == "patient_id") {
    set.seed(1)
    results <- results[sample(nrow(results)), ]
  }
  p <- ggplot(results, aes(UMAP_1, UMAP_2, colour = .data[[colour_by]])) +
    geom_point(size = 0.05, shape = 16) +
    facet_grid(facet_formula, labeller = label_both) +
    labs(title = title, colour = colour_by) +
    theme_bw() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), panel.grid = element_blank())
  if (colour_by == "cluster") {
    centers <- aggregate(results[, c("UMAP_1", "UMAP_2")],
                         by = results[, c("cluster", facet_vars)], FUN = median)
    p <- p + geom_text(data = centers, aes(label = cluster), colour = "black", size = 2.5) +
      NoLegend()
  } else {
    p <- p + guides(colour = guide_legend(override.aes = list(size = 3)))
  }
  p
}

##Harmony (cluster-level correction of the PCA embedding)
harmony_integration <-function(data,theta,lambda,max_iter,sigma){
   set.seed(42)
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
    
    #Return only the light table: one row per cell
    umap <- Embeddings(data, reduction = "umap.harmony")
    data.frame(cell       = rownames(umap),
               UMAP_1     = umap[, 1],
               UMAP_2     = umap[, 2],
               cluster    = data$clusters_harmony,
               patient_id = data$patient_id,
               params     = params_text,
               theta = theta, lambda = lambda, max_iter = max_iter, sigma = sigma,
               row.names  = NULL)
 }
     

#test_df <- harmony_integration(data, theta = 2, lambda = 1, max_iter = 10, sigma = 0.1)
#plot_umap_run(test_df, "umapharmony_")

#Each section has its own grid variable (harmony_grid / cca_grid), so running one section
#on its own can never pick up the other section's parameters.
harmony_grid <- expand.grid(theta = c(0,0.5, 1, 2, 4), lambda = 1, max_iter = 10, sigma = 0.1)
harmony_grid   # print it: 5 rows = 5 runs

harmony_results <- vector("list", nrow(harmony_grid))

for (i in 1:nrow(harmony_grid)) {
  harmony_results[[i]] <- harmony_integration(data, theta = harmony_grid$theta[i], lambda = harmony_grid$lambda[i],
                                              max_iter = harmony_grid$max_iter[i], sigma = harmony_grid$sigma[i])
  harmony_results[[i]]$run <- i
  gc()
}
harmony_results <- do.call(rbind, harmony_results)

#Save the table so the figures can be redone without re-running the integrations
dir.create("Module1 (Exploratory Analysis)/intermediate", showWarnings = FALSE)
qs_save(harmony_results, "Module1 (Exploratory Analysis)/intermediate/s1_harmony_grid_results.qs2")

#One PNG per run
for (i in 1:nrow(harmony_grid)) {
  umap_plot <- plot_umap_run(harmony_results[harmony_results$run == i, ], "umapharmony_")
  file_name <- file.path(harmony_dir, paste0("UMAP_theta", harmony_grid$theta[i],
                                             "_lambda", harmony_grid$lambda[i], "_max.iter.harmony", harmony_grid$max_iter[i],
                                             "_sigma", harmony_grid$sigma[i], ".png"))
  ggsave(filename = file_name, plot = umap_plot, width = 10, height = 5, dpi = 300)
}
rm(umap_plot)

#Combined overview: top row clusters, bottom row patients, one column per theta
harmony_overview <- patchwork::wrap_plots(
  plot_grid_overview(harmony_results, ~ theta, "cluster",    "Clusters"),
  plot_grid_overview(harmony_results, ~ theta, "patient_id", "Patients"),
  ncol = 1) +
  patchwork::plot_annotation(
    title    = "Harmony: efecto de theta",
    subtitle = paste0("Fijos: lambda = ", harmony_grid$lambda[1], " | sigma = ", harmony_grid$sigma[1], " | max_iter = ", harmony_grid$max_iter[1]))
ggsave(filename = "Module1 (Exploratory Analysis)/results/UMAP/harmony_grid_theta.png",
       plot = harmony_overview, width = 4 * nrow(harmony_grid), height = 9, dpi = 200)
rm(harmony_overview)
gc()

##Final Harmony run ON THE OBJECT. The grid above only kept light tables, so without this
##step data has no "harmony" reduction, no clusters_harmony and no umap.harmony (needed by s2).
##CHOOSE theta from harmony_grid_theta.png (2 is the Harmony default).
harmony_theta <- 2
set.seed(42)
data <- IntegrateLayers(object = data, method = HarmonyIntegration,
                        orig.reduction = "pca", new.reduction = "harmony",
                        normalization.method = "SCT", verbose = FALSE,
                        theta = harmony_theta, lambda = 1, max.iter.harmony = 10, sigma = 0.1)
data <- FindNeighbors(data, reduction = "harmony", dims = 1:15,
                      graph.name = c("harmony_nn", "harmony_snn"))
data <- FindClusters(data, graph.name = "harmony_snn", resolution = 0.5,
                     cluster.name = "clusters_harmony")
data <- RunUMAP(data, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony")


##CCA (Seurat anchors: cell-level mutual nearest neighbours in a shared CCA space)
##This is the method used in the original tumour analysis (Tumor_Analysis_Code.R).
##It is slower and uses more memory than Harmony, so it is done only ONCE here.
data <- IntegrateLayers(object = data, method = CCAIntegration,
                        orig.reduction = "pca", new.reduction = "integrated.cca",
                        normalization.method = "SCT", verbose = FALSE)

data <- FindNeighbors(data, reduction = "integrated.cca", dims = 1:15,
                      graph.name = c("cca_nn", "cca_snn"))
data <- FindClusters(data, graph.name = "cca_snn", resolution = 0.5,
                     cluster.name = "clusters_cca")

cca_integration <-function(data, n_neighbors, minimum_distance){
  set.seed(42)
  data <- RunUMAP(data, reduction = "integrated.cca", dims = 1:15, reduction.name = "umap.ccaintegration", n.neighbors =n_neighbors, min.dist = minimum_distance)
  
  #A text to write down the parameters values
  params_text <- paste0("minumum_distance = ",minimum_distance, " | num_neighbors = ", n_neighbors)
  
  #Return only the light table: one row per cell
  umap <- Embeddings(data, reduction = "umap.ccaintegration")
  data.frame(cell       = rownames(umap),
             UMAP_1     = umap[, 1],
             UMAP_2     = umap[, 2],
             cluster    = data$clusters_cca,
             patient_id = data$patient_id,
             params     = params_text,
             n_neighbors = n_neighbors, minimum_distance = minimum_distance,
             row.names  = NULL)
}


# test_plot <- cca_integration(data, n_neighbors = 30L, minimum_distance = 0.3 )
# test_plot


cca_grid <- expand.grid(n_neighbors = c(10L,20L,30L,40L,50L), minimum_distance =c(0.1,0.2,0.3,0.4,0.5))
cca_grid   # print it: 25 rows = 25 runs
cca_results <- vector("list", nrow(cca_grid))

for (i in 1:nrow(cca_grid)) {
  cca_results[[i]] <- cca_integration(data, n_neighbors = cca_grid$n_neighbors[i], minimum_distance = cca_grid$minimum_distance[i])
  cca_results[[i]]$run <- i
  gc()
}
cca_results <- do.call(rbind, cca_results)

#Save the table so the figures can be redone without re-running the UMAPs:
#cca_results <- qs_read(file.path(intermediate_dir, "s1_cca_grid_results.qs2"))
qs_save(cca_results, file.path(intermediate_dir, "s1_cca_grid_results.qs2"))

#One PNG per run
for (i in 1:nrow(cca_grid)) {
  umap_plot <- plot_umap_run(cca_results[cca_results$run == i, ], "umapccaintegration_")
  file_name <- file.path(cca_dir, paste0("UMAP_minimum_distance", cca_grid$minimum_distance[i],
                                         "_num_of_neighbors", cca_grid$n_neighbors[i], ".png"))
  ggsave(filename = file_name, plot = umap_plot, width = 10, height = 5, dpi = 300)
}
rm(umap_plot)

#Final CCA figure: two 5x5 grids instead of one 14 x 125 inch image.
#Rows = minimum_distance, columns = n_neighbors. The labels come from the table itself.
cca_grid_clusters <- plot_grid_overview(cca_results, minimum_distance ~ n_neighbors, "cluster",
                                        "CCA clusters: effect of number of neighbors and minimum distance")
cca_grid_patients <- plot_grid_overview(cca_results, minimum_distance ~ n_neighbors, "patient_id",
                                        "CCA patients: effect of number of neighbors and minimum distance")

ggsave(filename = "Module1 (Exploratory Analysis)/results/CCA/cca_grid_clusters.png",
       plot = cca_grid_clusters, width = 17.5, height = 17.5, dpi = 150)
ggsave(filename = "Module1 (Exploratory Analysis)/results/CCA/cca_grid_patients.png",
       plot = cca_grid_patients, width = 19, height = 17.5, dpi = 150)
rm(cca_grid_clusters, cca_grid_patients)
gc()

##Final CCA UMAP, with the setting chosen from the grid
data <- RunUMAP(data, reduction = "integrated.cca", dims = 1:15, reduction.name = "umap.cca",
                n.neighbors = 30L, min.dist = 0.3)


#Diagnostics: are the clusters biology or technical/patient effects?
##Marker genes and annotation are done in s2_script.R, on the object saved at the end.

##1. Cluster x patient: which fraction of each cluster comes from each patient
cl_pat <- table(cluster = data$clusters_cca, patient = data$patient_id)
cl_pat_frac <- round(prop.table(cl_pat, 1), 3)
cl_pat
cl_pat_frac
write.csv(as.data.frame.matrix(cl_pat),      file.path(diag_dir, "s1_cluster_x_patient_counts.csv"))
write.csv(as.data.frame.matrix(cl_pat_frac), file.path(diag_dir, "s1_cluster_x_patient_frac.csv"))

##Patient-specific clusters: one patient contributes more than 80% of the cells
names(which(apply(cl_pat_frac, 1, max) > 0.8))

pat_comp <- ggplot(as.data.frame(cl_pat), aes(cluster, Freq, fill = patient)) +
  geom_col(position = "fill") +
  labs(x = "CCA cluster", y = "Fraction of cells", fill = "patient_id") +
  theme_classic()
ggsave(file.path(diag_dir, "s1_cluster_patient_composition.png"), pat_comp, width = 8, height = 5, dpi = 300)

##2. QC by cluster. What to look for:
##   low-quality cells   -> low nFeature_RNA and/or high percent.mt (cluster 2?)
##   B/myeloid doublets  -> high B_frac / myeloid_frac (cluster 9?)
##   dissociation stress -> high dissoc_score1 (cluster 10?)
##   ambient RNA         -> high ambient_score1 spread over many cells
qc_cols <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "dissoc_score1", "ambient_score1", "dbl_score")
qc_by_cluster <- data@meta.data %>%
  group_by(clusters_cca) %>%
  summarise(n_cells = n(), across(all_of(qc_cols), median),
            dbl_frac     = mean(dbl_class == "doublet"),
            B_frac       = mean(B_doublet),
            myeloid_frac = mean(myeloid_doublet))
qc_by_cluster
write.csv(qc_by_cluster, file.path(diag_dir, "s1_qc_by_cluster.csv"), row.names = FALSE)

qc_vln <- VlnPlot(data, features = qc_cols, group.by = "clusters_cca", pt.size = 0, ncol = 2)
ggsave(file.path(diag_dir, "s1_qc_by_cluster.png"), qc_vln, width = 12, height = 10, dpi = 300)
rm(qc_vln)

##3. Integration metrics
##iLISI = effective number of patients in the neighbourhood of each cell (1 = a single patient,
##maximum = number of patients). Computed on the corrected embeddings (not on the UMAP).
##Careful: tumour cells are partly patient-specific for real (different clones), so the
##highest iLISI is not automatically the best one; too much mixing can mean over-correction.
emb <- c(PCA = "pca", Harmony = "harmony", CCA = "integrated.cca")
ilisi <- do.call(rbind, lapply(names(emb), function(m) {
  data.frame(method  = m,
             cluster = data$clusters_cca,
             iLISI   = lisi::compute_lisi(Embeddings(data, emb[[m]])[, 1:15],
                                          data@meta.data, "patient_id")$patient_id)
}))
ilisi$method <- factor(ilisi$method, levels = names(emb))
ilisi %>% group_by(method) %>% summarise(median_iLISI = median(iLISI))
ilisi_by_cluster <- ilisi %>% group_by(method, cluster) %>% summarise(median_iLISI = median(iLISI), .groups = "drop")
write.csv(ilisi_by_cluster, file.path(diag_dir, "s1_iLISI_by_cluster.csv"), row.names = FALSE)

ilisi_plot <- ggplot(ilisi, aes(method, iLISI)) +
  geom_boxplot(outlier.size = 0.1) +
  labs(x = NULL, y = "iLISI (patient_id)") +
  theme_classic()
ggsave(file.path(diag_dir, "s1_iLISI_methods.png"), ilisi_plot, width = 5, height = 5, dpi = 300)
rm(ilisi)

##ARI: how similar two clusterings are (1 = identical, 0 = as similar as by chance)
mclust::adjustedRandIndex(data$clusters_harmony, data$clusters_cca)

##ARI of every theta of the Harmony grid against CCA (uses the saved grid table, no re-run)
if (!exists("harmony_results")) harmony_results <- qs_read(file.path(intermediate_dir, "s1_harmony_grid_results.qs2"))
cca_cl <-setNames(as.character(data$clusters_cca), colnames(data))
ari_theta <- harmony_results %>%
  group_by(theta) %>%
  summarise(ARI_vs_CCA = mclust::adjustedRandIndex(cluster, cca_cl[cell]))
ari_theta
write.csv(ari_theta, file.path(diag_dir, "s1_ARI_harmony_theta_vs_cca.csv"), row.names = FALSE)

##Cluster correspondence Harmony vs CCA
write.csv(as.data.frame.matrix(table(Harmony = data$clusters_harmony, CCA = data$clusters_cca)),
          file.path(diag_dir, "s1_clusters_harmony_vs_cca.csv"))

##4. UMAP faceted by patient: all cells in grey, the cells of that patient coloured by cluster
u <- Embeddings(data, reduction = "umap.cca")
colnames(u) <- c("UMAP_1", "UMAP_2")
u <- cbind(as.data.frame(u), data@meta.data[, c("patient_id", "clusters_cca")])
umap_by_patient <- ggplot(u, aes(UMAP_1, UMAP_2)) +
  geom_point(data = u[, c("UMAP_1", "UMAP_2")], colour = "grey85", size = 0.05, shape = 16) +
  geom_point(aes(colour = clusters_cca), size = 0.1, shape = 16) +
  facet_wrap(~ patient_id) +
  labs(colour = "CCA cluster") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_classic()
ggsave(file.path(diag_dir, "s1_UMAP_cca_by_patient.png"), umap_by_patient, width = 14, height = 10, dpi = 200)
rm(u, umap_by_patient)
gc()


#Save the integrated object for s2_script.R
##Large file (several GB): it is in .gitignore, sync it to SurfDrive instead of git.
stopifnot(all(c("clusters_harmony", "clusters_cca") %in% colnames(data@meta.data)),
          all(c("harmony", "integrated.cca", "umap.harmony", "umap.cca") %in% Reductions(data)))
qs_save(data, file.path(intermediate_dir, "s1_integrated_harmony_cca.qs2"))
