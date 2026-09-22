set.seed(42)
#Download the necessary libraries
library(dplyr)
library(Seurat)
library(qs2)
library(ggplot2)

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
#We'll keep the cells with at leats 500 genes detected and maximal 7500(nFeatureRNA), and below a 10% percentage of mithocondrial genes. 

# NOTE (2026-09-22): subset() silently failed to trim the RNA assay's counts
# matrix (it updated meta.data correctly but left ncol(RNA) unchanged at 45,747)
# because the ADT assay held un-joined, per-sample split layers (counts.5/6/8/9/16-21)
# that broke Seurat's internal per-assay subsetting logic.
# Fix: ADT was unused, so it was dropped entirely (data[["ADT"]] <- NULL) before
# subsetting, restoring a single-layer RNA-only object on which subset() and
# NormalizeData() behave correctly.

data[["ADT"]] <- NULL

#Confirm that we have only RNA left
Assays(data)
data[["RNA"]]

data_filtered<- subset(data, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 10)

##Sanity check
ncol(data_filtered)
nrow(data_filtered@meta.data)
dim(data_filtered[["RNA"]])

#Normalization of the data
##From the tutorial (https://satijalab.org/seurat/articles/pbmc3k_tutorial),we follow the instructions and proccess the data like that

data_normalized<- NormalizeData(data_filtered)


