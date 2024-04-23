# Title: Create MDS plot
# Purpose: Visualise relationship of VSMC transcriptome to
# other publicly available RNASeq data
# Author: Charles Solomon
# Date: 21/02/2024
# Notes: (1)Run in R 4.1.0, (2)Save intermediate objects to disk to save time


# Load libraries ----------------------------------------------------------
library(data.table)
library(tidyverse)
library(DESeq2)
library(randomcoloR)
library(ggsci)
library(ggforce)

# Import dataset ----------------------------------------------------------
count <- data.frame(fread("mdsFiles/merged_expression_count_4_MDS_plot.txt"))
design <- data.frame(fread("mdsFiles/merged_expression_count_sample_attributes_4_MDS_plot.txt"))

rownames(count) <- count$GeneID
count$GeneID <- NULL
rownames(design) <- design$sampID
design$sampID <- NULL
design$alias <- factor(design$alias)


# all OK?
all(rownames(design) == colnames(count))

# Create a DESeq2 object --------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = design, design = ~ alias)

# filter genes whose sum of expression less than sum of samples
dds <- dds[ rowSums(counts(dds)) > ncol(dds), ]
save(dds, file = "mdsFiles/ddsObject.RData")

#load("dataFiles/ddsObject.RData")

# Variance transformation
vsd <- vst(dds, blind = FALSE)
save(vsd, file = "mdsFiles/vsdObject.RData")

load("mdsFiles/vsdObject.RData")

pcaData <- plotPCA(vsd, intgroup=c("alias"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=alias)) +
  geom_point(size=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_colour_manual(values=palette1)+
  theme_bw()

ggsave(paste0("mdsFiles/pcaPlot_", Sys.Date(), ".png"))



# Calculate euclidean distance
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
save(sampleDistMatrix, file = "mdsFiles/vsdDistMat.RData")

#load("dataFiles/vsdDistMat.RData")

rownames(sampleDistMatrix) <- vsd$alias
colnames(sampleDistMatrix) <- NULL

# Get MDS
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
save(mds, file = "mdsFiles/msdObject.RData")


###############################################################################
load("dataFiles/msdObject.RData")

# Plot MDS

palette <- c("#EC74B3", "#75E4D5", "#8B5085", "#7CAD71", "#E3D242", "#66E455", "#64E59E", "#7B6A6E", "#C79B42",
             "#E898E2", "#947BDB", "#D9E07E", "#DA72E6", "#62CAE0", "#E17648", "#A0B4AA", "#758BAE", "#E5A1C3", "#BBE4E8",
             "#E7E7D3", "#B3EC4F", "#C3E9BC", "#C2A6DF", "#94C4EF", "#7955DD", "#DA3CA3", "#9432E2", "#DDA6A5", "#6895E8",
             "#D9B989", "#DA5E74", "#AFE990", "#E03CE9")
palette1 <- sample(palette, 10)
palette1 <- c("#947BDB", "#EC74B3", "#B3EC4F", "#D9B989", "#E3D242", "#E898E2", "#DDA6A5", "#A0B4AA", "#7B6A6E", "#B3EC4F")

ggplot(mds, aes(x = `1`, y = `2`, color = alias)) +
  geom_mark_ellipse(aes(fill = alias, label = alias)) +
  geom_point(size = 1) + coord_fixed() +
  labs(title = "", x = "MDS Coordinate 1", y = "MDS Coordinate 2",
       color='Cells')+
  #scale_color_npg() +
  scale_colour_manual(values=palette1)+
  theme_bw()

ll <- design[!duplicated(design$alias), c(4,1,2)]
write.csv(ll, "mdsFiles/sampDetails.csv", row.names = F)


ggsave(paste0("mdsFiles/mdsPlot_", Sys.Date(), ".pdf"))
ggsave(paste0("mdsFiles/mdsPlot_", Sys.Date(), ".png"))
ggsave(paste0("mdsFiles/mdsPlot_", Sys.Date(), ".tiff"))

# ggsave("plots/mdsPlot_woHUVEC.pdf")
# ggsave("plots/mdsPlot_woHUVEC.png")

# Get variance explained
# The stress value in multidimensional scaling (MDS) is a measure of how well the distances in the lower-dimensional space represent the original dissimilarities. It is not directly convertible to variance explained in the way that, for example, R-squared is in linear regression.
# However, you can get a sense of the quality of the representation by considering the stress in relation to the maximum possible stress. The maximum stress occurs when the dissimilarities in the lower-dimensional space are completely unrelated to the original dissimilarities.
# The proportion of stress reduction can be calculated as follows:
# Proportion of Stress Reduction=Maximum Stress−Final StressMaximum Stress\text{Proportion of Stress Reduction} = \frac{\text{Maximum Stress} - \text{Final Stress}}{\text{Maximum Stress}}Proportion of Stress Reduction=Maximum StressMaximum Stress−Final Stress​
# In R, the metaMDS function in the vegan package returns a stress value. You can use this value to calculate the proportion of stress reduction:

# Install and load the 'vegan' package
# install.packages("vegan")
library(vegan)

# Perform non-metric MDS
mds_result <- metaMDS(sampleDistMatrix)

# Extract stress value
stress_value <- mds_result$stress

# Print stress value
print(paste("Stress value:", stress_value))
