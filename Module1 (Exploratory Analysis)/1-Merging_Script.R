# Here we will combine the Juan and Austin data 

#Read Data
##Nature data
data_nature<-qs_read("Inputs/External/stripped_s1_harmony_tumor_cells.qs2")
data

##Austin files
p45 <- readRDS("Inputs/External/Patients_Austin/s6745_S1.Rds")
p37 <- readRDS("Inputs/External/Patients_Austin/s6745_S2.Rds") 
p35 <- readRDS("Inputs/External/Patients_Austin/s6950_S1.Rds")
p44 <- readRDS("Inputs/External/Patients_Austin/s6950_S2.Rds")

##Combine the Austin files
sc_list <- list(s6745_s1 = p45, s6745_s2 = p37, s6950_s1 = p35, s6950_s1 = p44)

merged <- merge(
  x = sc_list[[1]],
  y = sc_list[-1],
  add.cell.ids = names(sc_list),   # avoid duplicated barcodes
  project = "NKI"
)

#We erase some columns from the datasets imported from Austin for merging purposes
keep <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "sample_id",
          "patient_id", "combi_id", "percent.mito", "barcode")

# 1- AUSTIN DATASET
pac_map <- c("6745_S1" = "P45",
             "6745_S2" = "P37",
             "6950_S1" = "P35",   
             "6950_S2" = "P44")

md <- merged@meta.data
md$patient_id   <- unname(pac_map[as.character(md$orig.ident)])
md$sample_id    <- ifelse(md$DF.classification == "Singlet", "SG_Singlets", NA)
md$percent.mito <- md$mito
md$combi_id     <- paste(md$patient_id, md$sample_id, sep = "_")  

bc <- regmatches(rownames(md), regexpr("[ACGT]+_[0-9]+$", rownames(md)))
md$barcode <- sub("_", "-", bc)

md$orig.ident <- "austin"
merged@meta.data <- md[, keep]

#Cell names with the same format as the original
merged <- RenameCells(merged, new.names = paste(merged$combi_id, merged$barcode, sep = "_"))

#2-Original dataset
data_nature$orig.ident <- "original"
data_nature@meta.data  <- data_nature@meta.data[, keep]

#3-Merge (just the RNA)
DefaultAssay(merged) <- "RNA"
DefaultAssay(data_nature) <- "RNA"
combined <- merge(DietSeurat(data_nature, assays = "RNA"),
                  DietSeurat(merged,      assays = "RNA"))

combined@meta.data <- combined@meta.data[, keep]


#4 Comprobations
colnames(combined@meta.data)
table(combined$orig.ident)
table(combined$patient_id, combined$orig.ident)
table(combined$sample_id,  combined$orig.ident)
any(duplicated(colnames(combined)))
summary(data_nature$percent.mito); summary(merged$percent.mito)

# 5 store
saveRDS(combined, file = "Inputs/External/combined.rds")

#Save session information
Info_script <-sessionInfo()

#Save environment
save.image("Module1 (Exploratory Analysis)/intermediate/merging.rdata")
