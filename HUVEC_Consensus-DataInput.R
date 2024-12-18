# 
# Charles Solomon
# 01/11/2022


# Load libraries and make settings-----------------------------------------
library(data.table)
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Import gene expression data for tnfa and untr -------------------------------
tnfaData <-  data.frame(fread("exprFiles/tnfa_pilot_BRBseq_geneExpr_qtlReady.txt"))
untrData <-  data.frame(fread("exprFiles/untreated_pilot_BRBseq_geneExpr_qtlReady.txt"))

# Select only genes present in both datasets
allData <- merge(tnfaData, untrData, by = "geneID")
tnfaData <- tnfaData[tnfaData$geneID %in% allData$geneID, ]
untrData <- untrData[untrData$geneID %in% allData$geneID, ]


# We work with two sets:
nSets = 2;
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("TNFa HUVEC", "UNTR HUVEC")
shortLabels = c("Tnfa", "Untr")
# Form multi-set expression data
multiExpr = vector(mode = "list", length = nSets)

multiExpr[[1]] = list(data = as.data.frame(t(tnfaData[, -c(1)])));
names(multiExpr[[1]]$data) = tnfaData$geneID;
rownames(multiExpr[[1]]$data) = names(tnfaData)[-c(1)];
multiExpr[[2]] = list(data = as.data.frame(t(untrData[-c(1)])));
names(multiExpr[[2]]$data) = untrData$geneID;
rownames(multiExpr[[2]]$data) = names(untrData)[-c(1)];
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)


# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK


if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}


# Plot Sample Tree
sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

pdf(file = "wgcnaPlots/SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
dev.off();

# # Cut Tree if necessary
# 
# # Choose the "base" cut height for the tnfa data set
# baseHeight = 157
# # Adjust the cut height for the untr data set for the number of samples
# cutHeights = c(157, 157);
# # Re-plot the dendrograms including the cut lines
# pdf(file = "Plots/SampleClustering.pdf", width = 12, height = 12);
# par(mfrow=c(2,1))
# par(mar = c(0, 4, 2, 0))
# for (set in 1:nSets)
# {
#   plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
#        xlab="", sub="", cex = 0.7);
#   abline(h=cutHeights[set], col = "red");
# }
# dev.off();
# 
# 
# for (set in 1:nSets)
# {
#   # Find clusters cut by the line
#   labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set])
#   # Keep the largest one (labeled by the number 1)
#   keep = (labels==1)
#   multiExpr[[set]]$data = multiExpr[[set]]$data[keep, ]
# }
# collectGarbage();
# Check the size of the leftover data
exprSize = checkSets(multiExpr)
exprSize

# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;


# Save input data
save(multiExpr, nGenes, nSamples, setLabels, shortLabels, exprSize, 
     file = "wgcnaFiles/HUVEC_Consensus-dataInput.RData")


