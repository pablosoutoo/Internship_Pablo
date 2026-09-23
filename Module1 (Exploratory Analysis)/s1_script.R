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

# NOTE (2026-09-23): subset() is unreliable on this object. It trims meta.data
# correctly (37,407 QC-passing cells) but leaves the RNA assay's counts matrix at
# the original 45,747 cells, so NormalizeData() fails with a dimension mismatch.
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

#For the write-up: record the versions and assay class involved in the subset() bug
class(data[["RNA"]])
packageVersion("Seurat")
packageVersion("SeuratObject")

##QC-passing barcodes, taken from meta.data
qc_pass <- with(data@meta.data, nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 10)
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
data_filtered <- CreateSeuratObject(counts = counts_filtered, meta.data = meta_filtered,
                                    assay = "RNA", min.cells = 0, min.features = 0)
rm(counts, counts_filtered, meta_filtered)

##Sanity check: all three must agree (expected 28,126 genes x 37,407 cells)
ncol(data_filtered)
nrow(data_filtered@meta.data)
dim(data_filtered[["RNA"]])
stopifnot(
  ncol(data_filtered) == length(keep_cells),
  nrow(data_filtered@meta.data) == length(keep_cells),
  ncol(data_filtered[["RNA"]]) == length(keep_cells),
  identical(colnames(data_filtered), rownames(data_filtered@meta.data))
)

#Normalization of the data
##From the tutorial (https://satijalab.org/seurat/articles/pbmc3k_tutorial),we follow the instructions and proccess the data like that

data_normalized<- NormalizeData(data_filtered)


