
# This script takes one argument -  the prefix name of the plink pca output
aha <- commandArgs(trailingOnly = TRUE)

aha2 <- paste0("genoFiles/", aha, ".eigenvec")
aha3 <- paste0("covFiles/", aha, "_20PCs.txt")

# read in the eigenvectors, produced in PLINK
eigenvec <- read.table(aha2, header = FALSE, skip=0, sep = ' ')
rownames(eigenvec) <- eigenvec[,2]
eigenvec <- eigenvec[,3:ncol(eigenvec)]
colnames(eigenvec) <- paste("PC", c(1:20), sep = '')
eigenvec <- data.frame(t(eigenvec))
eigenvec <- data.frame(Covariates = rownames(eigenvec), eigenvec)

write.table(eigenvec, aha3, sep = "\t")

