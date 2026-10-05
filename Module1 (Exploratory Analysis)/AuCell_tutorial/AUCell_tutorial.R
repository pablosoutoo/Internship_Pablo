#Install package
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
# To support paralell execution:
BiocManager::install(c("doMC", "doRNG","doSNOW"))
# For the main example:
BiocManager::install(c("mixtools", "SummarizedExperiment"))
# For the examples in the follow-up section of the tutorial:
BiocManager::install(c("DT", "plotly", "NMF", "d3heatmap",
                       "dynamicTreeCut", "R2HTML", "Rtsne", "zoo"))

#Expression matrix
##Download the data
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
library(GEOquery)

geoFile <- "GSE60361_C1-3005-Expression.txt.gz"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60361/suppl/GSE60361_C1-3005-Expression.txt.gz", destfile = geoFile)

library(data.table)
exprMatrix <-fread(geoFile, sep="\t")
geneNames <- unname(unlist(exprMatrix[,1, with=FALSE]))
exprMatrix <- as.matrix(exprMatrix[,-1, with=FALSE])
rownames(exprMatrix) <- geneNames
exprMatrix <- exprMatrix[unique(rownames(exprMatrix)),]
dim(exprMatrix)
exprMatrix[1:5,1:4]

##Remove file(s) downloaded
file.remove(geoFile)

#Convert to sparse
library(Matrix)
exprMatrix <- as(exprMatrix, "dgCMatrix")

#Save for future use
mouseBrainExprMatrix <- exprMatrix
save(mouseBrainExprMatrix, file="exprMatrix_MouseBrain.RData")

##To speed-up the execution of the tutorial, we will use only 5000 random genes from this dataset.

# load("exprMatrix_MouseBrain.RData")
set.seed(333)
exprMatrix <- mouseBrainExprMatrix[sample(rownames(mouseBrainExprMatrix), 5000),]


#Genesets - The gene sets or signatures to test on the cells
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSEABase")
BiocManager::install("AUCell")
library(AUCell)
library(GSEABase)

gmtFile <- paste(file.path(system.file('examples', package='AUCell')), "geneSignatures.gmt", sep="/")
geneSets <- getGmt(gmtFile)
geneSets

geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix)) 
cbind(nGenes(geneSets))

geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets), " (", nGenes(geneSets) ,"g)", sep=""))


#For the example, we also add a few sets of random genes and 100 genes expressed in many vells
##Random
set.seed(321)
extraGeneSets <-c(
  GeneSet(sample(rownames(exprMatrix), 50), setName = "Random (50g)"),
  GeneSet(sample(rownames(exprMatrix), 500), setName ="Random (500g)"))

countsPerGene <- apply(exprMatrix, 1, function(x) sum(x>0))
# Housekeeping-like
extraGeneSets <- c(extraGeneSets,
                   GeneSet(sample(names(countsPerGene)[which(countsPerGene>quantile(countsPerGene, probs=.95))], 100), setName="HK-like (100g)"))

geneSets <- GeneSetCollection(c(geneSets,extraGeneSets))
names(geneSets)

#SCORE GENE SIGNATUES
#The function AUCell_run calculates the signature enrichment scores on each cell

cells_AUC <- AUCell_run(exprMatrix, geneSets)
save(cells_AUC, file="cells_AUC.RData")

# The AUC scores are calculated in two steps: 1. Building the rankings (AUCell_buildRankings) and 2. Calculate enrichment (AUCell_calcAUC).
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=TRUE)
cells_rankings

save(cells_rankings, file="cells_rankings.RData")

#Calculate enrinchment for the gene signatures (AUC)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
save(cells_AUC, file="cells_AUC.RData")

#Exploration of distributions
set.seed(333)
par(mfrow=c(3,3)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 


warningMsg <- sapply(cells_assignment, function(x) x$aucThr$comment)
warningMsg[which(warningMsg!="")]

cells_assignment$Oligodendrocyte_Cahoy$aucThr$thresholds

#To obtain the threshold selected automatically for a given gene set(e.g. Oligodendrocytes)
cells_assignment$Oligodendrocyte_Cahoy$aucThr$selected

#Cells asigned at this threshold
oligodencrocytesAssigned <- cells_assignment$Oligodendrocyte_Cahoy$assignment
length(oligodencrocytesAssigned)
head(oligodencrocytesAssigned)

#Plotting the AUC histogram of a specific gene set, and setting a new threshold:
geneSetName <- rownames(cells_AUC)[grep("Oligodendrocyte_Cahoy", rownames(cells_AUC))]
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=0.25)
abline(v=0.25)

#Assigning cells to this new threshold
newSelectedCells <- names(which(getAUC(cells_AUC)[geneSetName,]>0.08))
length(newSelectedCells)
head(newSelectedCells)
