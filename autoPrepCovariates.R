# Prepare covariates
# Charles Solomon
# 15/07/2022

# This script takes three arguments 
#  - the prefix name of the plink pca output
#  - number of genotype PC to use
#  - number of peer factors to use

# Assign command line arguments
aha <- commandArgs(trailingOnly = TRUE)
aha1 <- aha[1]
aha2 <- as.numeric(aha[2])
aha3 <- as.numeric(aha[3])

sex <- read.table(paste0("covFiles/", aha1, "_Sex.txt"))
sex <- sex[, order(colnames(sex))]

pca <- read.table(paste0("covFiles/", aha1, "_20PCs.txt"))
pca <- pca[, order(colnames(pca))]

peer <- read.table(paste0("covFiles/", aha1, "_PeerFactor.txt"), header = T)
peer <- peer[, order(colnames(peer))]

cov <- rbind(sex, pca[1:aha2, ], peer[1:aha3, ])

write.table(cov, paste0("covFiles/", aha1, "_", aha2, "pc_", aha3,"peer.txt"), sep = "\t", row.names = F)

cat("Covariates file Ready  \n")
