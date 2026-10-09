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
violin_plot
ggsave(filename="Module1 (Exploratory Analysis)/results/s1_mitocond_pres.png", plot=violin_plot,width = 8, height = 6, dpi = 300)

##FeatureScatter <- Tipically used to visualize feature-feature relationship

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")

plot1
plot2

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


#Drop some genes that we are not interested in
drop <- c(grep("^(MT-|                      #Mitocondrial genes
               RP[SL]\\d|                   #Ribosomal proteins
               TR[ABDG][VJC]|               #V,J and C segmenst of the TCR (clonotipe specific and they would split the cells by clone)
               IG[HKL][VJC])",              #V,J and C segments of immunoglobulines
               rownames(data), value = TRUE),
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



#Elbow plot to see the influence of the number of dimensions in the reduction of standar deviation.
elbow_plot<-ElbowPlot(data)
elbow_plot
ggsave(filename="Module1 (Exploratory Analysis)/results/PCA/Elbow_plot.png", plot=elbow_plot,width = 8, height = 6, dpi = 300)


#Integration
#We integrate the same object with two different methods (Harmony and CCA) so we can compare them afterwards.
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
  
  #Where to put the tag of each cluster. Takes the median of the XY coordinates of all it's cells.
  centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ cluster, data = df, FUN = median)
  
  #Umap per clusters
  plot_clusters <- ggplot(df, aes(UMAP_1, UMAP_2, colour = cluster)) +
    geom_point(size = pt_size, shape = 16) +
    geom_text(data = centers, aes(label = cluster), colour = "black", size = 4) +
    labs(x = x_lab, y = y_lab) +
    theme_classic() +
    NoLegend() 
  
  #Umap per patients (to see if they are equally distributed)
  set.seed(42)
  df_shuffled <- df[sample(nrow(df)), ] #This is done so the points are painted on a random order and a patient does not overlap another
  plot_patients <- ggplot(df_shuffled, aes(UMAP_1, UMAP_2, colour = patient_id)) +
    geom_point(size = pt_size, shape = 16) +
    labs(x = x_lab, y = y_lab, colour = "patient_id") +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    theme_classic() 
  
  #Combine both panels in one figure
  patchwork::wrap_plots(plot_clusters, plot_patients) +
    patchwork::plot_annotation(title = df$params[1],
                               theme = theme(plot.title = element_text(size = 11)))
}

##Overview of a whole grid in ONE ggplot, faceted by the parameters (built from the table,
##so it cannot pick up a stale grid or plot_list). colour_by = "cluster" or "patient_id".
plot_grid_overview <- function(results, facet_formula, colour_by, title) {
  facet_vars <- all.vars(facet_formula) #Return a character vector containing all the names which occur in an expression or call.
  if (colour_by == "patient_id") {
    set.seed(42)
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
    ## Compute the k nearest neighbors for the data. The dimensions were chose based on the elbow plot.
    data <- FindNeighbors(data, reduction = "harmony", dims = 1:15,
                          graph.name = c("harmony_nn", "harmony_snn"))
    
    #Identification of clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm.
    data <- FindClusters(data, graph.name = "harmony_snn", resolution = 0.5,
                         cluster.name = "clusters_harmony")
    
    #Runs the Uniform Manifold Approximated and Projection (UMAP) dimensional reduction technique
    data <- RunUMAP(data, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony")
    
    #A text to write down the parameters values
    params_text <- paste0("theta = ",theta, " | lambda = ", lambda, " | max.iter.harmony =", max_iter, " | sigma = ", sigma)
    
    #Return only the light table: one row per cell. 
    #This was done because otherwise it used too much memory.
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
     
#Test done to see if the function was working

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
    title    = "Harmony: theta effect",
    subtitle = paste0("Fixed: lambda = ", harmony_grid$lambda[1], " | sigma = ", harmony_grid$sigma[1], " | max_iter = ", harmony_grid$max_iter[1]))
ggsave(filename = "Module1 (Exploratory Analysis)/results/UMAP/harmony_grid_theta.png",
       plot = harmony_overview, width = 4 * nrow(harmony_grid), height = 9, dpi = 200)
rm(harmony_overview)
gc()


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


#Save the data so it is not necessary to re-run the previous code
qs_save(data, "Module1 (Exploratory Analysis)/intermediate/s1_data_cca.qs2")


#Build the function of integration with CCA
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
             n_neighbors = n_neighbors, 
             minimum_distance = minimum_distance,
             row.names  = NULL)
}



# test_plot <- cca_integration(data, n_neighbors = 30L, minimum_distance = 0.3 )
# test_plot

#Build the grid to combine the different number of neighbors and minumum distance in just one plot.
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
cca_results <- qs_read(file.path(intermediate_dir, "s1_cca_grid_results.qs2"))
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



#Graphic to see the composition per patient in each cluster
comp <- data@meta.data %>%
  count(cluster = as.character(clusters_cca), patient = patient_id) %>%
  bind_rows(data@meta.data %>% count(patient = patient_id) %>% mutate(cluster = "All")) %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n), total = sum(n)) %>%
  ungroup() %>%
  mutate(cluster = factor(cluster, levels = c("All", levels(data$clusters_cca))))


p_comp <- ggplot(comp, aes(x = cluster, y = prop, fill = patient)) +
  geom_col(width = 0.8) +
  geom_text(data = distinct(comp, cluster, total),
            aes(x = cluster, y = 1.01, label = total),
            inherit.aes = FALSE, size = 3, vjust = 0) +          # nº de células encima de cada barra
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.06))) +
  labs(x = "CCA cluster", y = "Cells (%)", fill = "Patient",
       title = "Patient contribution to each CCA cluster") +
  theme_classic()


p_comp
ggsave("Module1 (Exploratory Analysis)/results/Annotation/composition_cluster_per_patient.png", plot = p_comp,width = 19, height = 17.5, dpi = 150 )

#FindAllMarkers
##SCTransform was run per patient, so the SCT assay holds one model per patient.
##FindMarkers refuses that until PrepSCTFindMarkers() puts all cells on the same scale.


data <- PrepSCTFindMarkers(data, assay = "SCT") #Given a merged object with multiple SCT models (in this case data)
                                                #this function uses minumum of the median UMI (calculated using the raw UMI counts) of 
                                                #individual objects to reverse the individual SCT regression model using minumum
                                                #of median UMI as the sequencing depth covariate.
Idents(data) <- "clusters_cca"
data.markers <- FindAllMarkers(data, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.2)
stopifnot(nrow(data.markers) > 0)   # if this stops, run warnings() to see why each cluster failed

data.markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

tab <- table(data$clusters_cca, data$patient_id)
round(prop.table(tab, 1), 2)          # composición de cada cluster
round(prop.table(table(data$patient_id)), 2)   # composición global, para comparar
tab

round(prop.table(tab, 2), 2)                                   # % de cada paciente en cada cluster
round(sweep(prop.table(tab, 1), 2, prop.table(table(data$patient_id)), "/"), 2)  # enriquecimiento

#Cell-population annotation of the CCA clusters
##Same idea as the Seurat PBMC tutorial (one name per cluster), but these are melanoma tumour
##cells, so the "populations" are tumour states (melanocytic, mitotic, neural crest-like,
##mesenchymal-like, IFN response, stress...; Pozniak et al. 2024) rather than immune cell types.
##The label is stored as a metadata column (cell_state_cca), so clusters_cca stays untouched.
annot_dir <- "Module1 (Exploratory Analysis)/results/Annotation"
dir.create(annot_dir, showWarnings = FALSE, recursive = TRUE)

##1. Evidence: top markers per cluster + canonical genes of each state
top_markers <- data.markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  mutate(delta_pct = pct.1 - pct.2) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC * delta_pct, n = 10, with_ties = FALSE) %>%
  ungroup()

canonical <- c("MKI67", "TOP2A",                          #cycling
               "MITF", "PMEL", "DCT", "MLANA", "TYR",     #melanocytic
               "SOX10", "NGFR", "AXL", "SOX9",            #neural crest-like / dedifferentiated
               "VIM", "SERPINE1", "FN1",                  #mesenchymal
               "VEGFA", "CA9", "NDUFA4L2",                #hypoxia
               "CXCL10", "STAT1", "B2M", "HLA-A",         #IFN response / MHC-I
               "CD74", "HLA-DRA",                         #MHC-II
               "HSPA6", "HSPA1A", "CDKN1A", "GDF15")      #stress (heat shock, p53)


Tumor_signatures<-read.csv("Inputs/External/Tumor_Signatures.csv", header=TRUE, row.names = 1)
 

gene_states <- Tumor_signatures %>%
  distinct(gene, Library, Cell_state) %>%
  group_by(gene) %>%
  summarise(Cell_state = paste(paste0(Cell_state, " (", Library, ")"), collapse = "; "),
            .groups = "drop")

top_markers_annot <- top_markers %>%
  ungroup() %>%
  left_join(gene_states, by = "gene")    

write.csv(top_markers_annot, file.path(annot_dir, "s1_top10_markers_cca.csv"), row.names = FALSE)


canonical <- intersect(canonical, rownames(data))
dot_cca <- DotPlot(data, features = canonical, group.by = "clusters_cca", assay = "SCT") +
  RotatedAxis() + ggtitle("CCA clusters - canonical genes")
dot_cca
ggsave(file.path(annot_dir, "s1_canonical_dotplot_cca.png"), dot_cca, width = 12, height = 6, dpi = 300)

##2. One name per cluster. Fill these in after looking at the dot plot and the marker table;
##clusters you are not sure about stay "Unassigned". Two clusters may share the same name.
AuCell_labels <- setNames(rep("Unassigned", nlevels(data$clusters_cca)), levels(data$clusters_cca))
AuCell_labels["0"]  <- "Neural crest-like"            #SOX10 high; NCMAP, SCN7A, SCRG1, ANGPTL7
AuCell_labels["1"]  <- "Melanocytic"                  #PMEL, MLANA, MITF, DCT, TYR
AuCell_labels["2"]  <- "Stress (ATF4 / amino acid)"   #ASNS, TRIB3, GDF15, CDKN1A (weak markers)
AuCell_labels["3"]  <- "IFN response"                 #GBP1/4, IFIT2, IFI44L, STAT1, B2M, HLA-A
AuCell_labels["4"]  <- "Mesenchymal-like (invasive)"  #MMP1, MMP3, IL11, SERPINB2, SERPINE1, INHBA
AuCell_labels["5"]  <- "Stress (hypoxia)"             #NDUFA4L2, VEGFA, CA9, MT3
AuCell_labels["6"]  <- "Melanocytic (pigmentation)"   #TYR, MITF, DCT high; low MHC-I
AuCell_labels["7"]  <- "Mesenchymal-like (TGFb/YAP)"  #FN1, TAGLN, CCN1, CCN2, DKK1
AuCell_labels["8"]  <- "Mitotic (G1/S)"               #E2F2, RRM2, MCM10, CDC45, CLSPN
AuCell_labels["9"]  <- "Antigen presentation (MHC-II)" #CD74, HLA-DRA; B-cell genes in ~5% of cells
AuCell_labels["10"] <- "Inflammatory (NF-kB)"         #CXCL10/11, CCL2, CXCL2, SELE, HSPA6
AuCell_labels["11"] <- "Mitotic (G2/M)"               #PLK1, CDC20, KIF20A, MKI67, TOP2A
AuCell_labels

data$cell_state_cca <- unname(AuCell_labels[as.character(data$clusters_cca)])
table(data$clusters_cca, data$cell_state_cca)

##3. Labelled UMAP (one final CCA UMAP, with the setting you chose from the grid)
data <- RunUMAP(data, reduction = "integrated.cca", dims = 1:15, reduction.name = "umap.cca",
                n.neighbors = 30L, min.dist = 0.3)
umap_states <- DimPlot(data, reduction = "umap.cca", group.by = "cell_state_cca",
                       label = TRUE, repel = TRUE, pt.size = 0.5) + NoLegend() +
  ggtitle("CCA clusters - cell states")
umap_states
ggsave(file.path(annot_dir, "s1_UMAP_cell_states_cca.png"), umap_states, width = 8, height = 6, dpi = 300)





#AUCELL
library(GEOquery)
library(AUCell)
library(data.table)

rna <- JoinLayers(data[["RNA"]])
exprMatrix <- LayerData(rna, layer = "counts")

#Gene Sets
##Before using the CSV, two filters
##Karras and Pozniak signatures are identical
##Remove Patient_specific A/B and Mitochondrial(low quality) -> They are not biological states
sig <- read.csv("Inputs/External/Tumor_Signatures.csv", row.names = 1) |>
  dplyr::filter(Library != "Pozniak", !grepl("Patient_specific|low_quality", Cell_state))
gs_list <- split(sig$gene, paste(sig$Library, sig$Cell_state))
gs_list <- lapply(gs_list, function(g) intersect(unique(g), rownames(exprMatrix)))
gs_list <- gs_list[lengths(gs_list) >= 10]


#Rankings and AUC
set.seed(42)
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats = TRUE)
cells_AUC <- AUCell_calcAUC(gs_list, cells_rankings,
                            aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))
qs_save(cells_AUC, file.path(intermediate_dir, "s1_cells_AUC.qs2"))


#Take a look at the results
auc_mat <- t(getAUC(cells_AUC))
colnames(auc_mat) <- make.names(colnames(auc_mat))
data <- AddMetaData(data, as.data.frame(auc_mat))

mean_auc <- aggregate(as.data.frame(auc_mat), list(cluster = data$clusters_cca), mean)
rownames(mean_auc) <- mean_auc$cluster
heatmap <- pheatmap::pheatmap(scale(as.matrix(mean_auc[, -1])))   # clusters × signatures
ggsave("Module1 (Exploratory Analysis)/results/Annotation/heatmap_signs_clus.png", plot = heatmap, width = 10, height = 5, dpi = 300)

#UMAP after AUCell
## One name per cluster. Fill these in after looking at the dot plot and the marker table;
##clusters you are not sure about stay "Unassigned". Two clusters may share the same name.
AuCell_labels <- setNames(rep("Unassigned", nlevels(data$clusters_cca)), levels(data$clusters_cca))
AuCell_labels["0"]  <- "Unassigned_1"            
AuCell_labels["1"]  <- "Melanocytic"                  
AuCell_labels["2"]  <- "Stress (ATF4 / amino acid)"   
AuCell_labels["3"]  <- "Immune/Antigen-presenting"                 
AuCell_labels["4"]  <- "Neural-crest-like/Undifferenciated"  
AuCell_labels["5"]  <- "Stress (hypoxia)"             
AuCell_labels["6"]  <- "Stressed/MSC-like"   
AuCell_labels["7"]  <- "Unassigned_2"  
AuCell_labels["8"]  <- "Mitosis"               
AuCell_labels["9"]  <- "Patient specific" 
AuCell_labels["10"] <- "IFN-responsive/invasive"         
AuCell_labels["11"] <- "Proliferating"               
AuCell_labels

data$AUCell_states <- unname(AuCell_labels[as.character(data$clusters_cca)])
table(data$clusters_cca, data$cell_state_cca)

##3. Labelled UMAP (one final CCA UMAP, with the setting you chose from the grid)
data <- RunUMAP(data, reduction = "integrated.cca", dims = 1:15, reduction.name = "umap.cca",
                n.neighbors = 30L, min.dist = 0.3)
umap_states <- DimPlot(data, reduction = "umap.cca", group.by = "clusters_cca",
                       label = TRUE, repel = TRUE, pt.size = 0.5) + NoLegend() +
  ggtitle("AuCell - cell states")
umap_states
ggsave(file.path(annot_dir, "s1_UMAP_cell_states_cca_AUCell.png"), umap_states, width = 8, height = 6, dpi = 300)


#See if the how the clusters are organized in melanocityc/immunitary

z <- scale(as.matrix(mean_auc[, -1]))   # clusters x signatures, z across clusters
states <- list(
  Melanocytic = c("Wouters.Melanocytic.cell.state", "Tsoi.Melanocytic",
                  "Rambow.MITFtargets", "Rambow.pigmentation"),
  Immune      = c("Karras.Antigen_presentation", "Rambow.Immune"))
sigs <- unlist(states)
row_state <- rep(names(states), lengths(states))

ComplexHeatmap::Heatmap(t(z[, sigs]), name = "z (mean AUC)",
                        row_split = factor(row_state, levels = names(states)),
                        cluster_rows = FALSE, cluster_row_slices = FALSE,
                        col = circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")))

state_score <- sapply(states, function(s) rowMeans(z[, s, drop = FALSE]))  # summary per cluster

#Most enriched signature per cluster
top1 <- data.frame(cluster   = rownames(z),
                   signature = colnames(z)[apply(z, 1, which.max)],
                   z         = round(apply(z, 1, max), 2))
top1

#Create a boxplot to see the distrubution of the signatures in each cluster
for (s in unique(top1$signature)) {
  top_clusters <- top1$cluster[top1$signature == s]          # cluster(s) where s is the most enriched
  df <- data.frame(auc     = auc_mat[colnames(data), s],     # same cell order as data
                   cluster = data$clusters_cca,
                   patient = data$patient_id)
  
  p <- ggplot(df, aes(cluster, auc, fill = patient)) +
    geom_boxplot(outlier.size = 0.3) +
    theme_bw() +
    labs(y = paste("AUC:", s), x = "Cluster",
         title = s, subtitle = paste("Top signature of cluster(s):", paste(top_clusters, collapse = ", ")))
  
  print(p)                                                    # para verlo en RStudio
  ggsave(file.path("Module1 (Exploratory Analysis)/results/Annotation/boxplots.png", paste0("box_", s, ".png")), p, width = 10, height = 4.5, dpi = 300)
}

#Save session information
Info_script <-sessionInfo()

#Save environment
save.image("Module1 (Exploratory Analysis)/intermediate/environment.rdata")
