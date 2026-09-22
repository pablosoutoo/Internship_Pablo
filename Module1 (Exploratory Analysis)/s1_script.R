set.seed()
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
volcano_plot<-VlnPlot(data, features =c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0, ncol=3)

ggsave(filename="Module1 (Exploratory Analysis)/results/s1_mitocond_pres.png", plot=volcano_plot,width = 8, height = 6, dpi = 300)
