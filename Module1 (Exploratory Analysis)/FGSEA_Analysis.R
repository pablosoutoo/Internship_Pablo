##FGSEA ANALYSIS

#Libraries
library(fgsea)
library(data.table)
library(ggplot2)
library(msigdbr)
library(qs2)
library(AUCell)

#Data
##Downlaad all the C2:CP collection for Human (Homo Sapiens)
c2_cp_human <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP")
data<-qs_read("Module1 (Exploratory Analysis)/intermediate/data_clusterized.qs2")

