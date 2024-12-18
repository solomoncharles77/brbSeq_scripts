# 
# Charles Solomon
# 23/08/2024


# Load libraries and make settings-----------------------------------------
library(data.table)
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Import gene expression data for untr and tnfa -------------------------------
tnfaData <-  data.frame(fread("exprFiles/tnfa_pilot_BRBseq_geneExpr_qtlReady.txt"))
untrData <-  data.frame(fread("exprFiles/untreated_pilot_BRBseq_geneExpr_qtlReady.txt"))

# Select only genes present in both datasets
allData <- merge(tnfaData, untrData, by = "geneID")
tnfaData <- tnfaData[tnfaData$geneID %in% allData$geneID, ]
untrData <- untrData[untrData$geneID %in% allData$geneID, ]


# Import gene expression and phenotype data -------------------------------
allExpr <- untrData

# Take a quick look at what is in the data set:
dim(allExpr)

# Convert expression data to WGCNA matrix ---------------------------------
datExpr0 <-  as.data.frame(t(allExpr[, -c(1)]))
colnames(datExpr0) <-  allExpr$geneID
rownames(datExpr0) <-  colnames(allExpr)[-c(1)]
gc()

# Check for excessive missing values --------------------------------------
gsg <-  goodSamplesGenes(datExpr0, verbose = 5)
gsg$allOK


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Cluster samples to detect possible outliers -----------------------------
sampleTree = hclust(dist(datExpr0), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


# # To cut or not to cut? ---------------------------------------------------
# # some samples seem like outliers so we attempted to cut them off
# abline(h = 26, col = "red");
# # Determine cluster under the line
# clust <-  cutreeStatic(sampleTree, cutHeight = 26, minSize = 1)
# table(clust)
# # clust 1 contains the samples we want to keep.
# keepSamples = (clust==1)
# datExpr = datExpr0[keepSamples, ]
### --------------------------------------------------------------------------

datExpr = datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# save plot
sizeGrWindow(12,9)
png(file = "wgcnaPlots/huvecUNTR_WGCNA_sampleClustering.png", width = 480, height = 480, units = "px");
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()


collectGarbage();


# Choose a set of soft-thresholding powers --------------------------------
powers <-  c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft <-  pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
png(file = "wgcnaPlots/huvecUNTR_WGCNA_softThresholdPower.png", width = 480, height = 380)
#sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.99,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

# based on the plots we chose the soft thresholding power 3

save(datExpr, sft, sampleTree,  file = "wgcnaFiles/huvecUNTR_WGCNA_dataInput.RData" )

#############################################################