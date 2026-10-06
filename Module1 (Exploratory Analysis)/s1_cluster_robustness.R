set.seed(42)

#This Script will be used as an exploratory analysis of the cluster distribution 

#Output folders
harmony_dir      <- "Module1 (Exploratory Analysis)/results/UMAP/harmony_grid"
cca_dir          <- "Module1 (Exploratory Analysis)/results/CCA/cca_grid"
intermediate_dir <- "Module1 (Exploratory Analysis)/intermediate"

#Libraries
library(Seurat)
library(qs2)
library(clustree)
library(ggplot2)
library(mclust)

#Data
data <- qs_read("Module1 (Exploratory Analysis)/intermediate/s1_data_cca.qs2")


#Resolution iterations
data <- FindClusters(data, graph.name = "cca_snn",
                     resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2))

cluster_tree<-clustree(data, prefix = "cca_snn_res.")
ggsave(filename="Module1 (Exploratory Analysis)/intermediate/cluster_tree.png", plot=cluster_tree,width = 8, height = 6, dpi = 300)


