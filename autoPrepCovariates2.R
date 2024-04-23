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


pca <- read.table(paste0("huvecImputeCov/", aha1, "_20PCs.txt"))
colnames(pca) <- gsub(".*S", "E", colnames(pca))
pca <- pca[, order(colnames(pca))]

sex <- read.delim("../../../shared/HUVEC_genotype/Imputed_Genotypes_Data_All_HUVEC_Samples2.fam", header = F)
sex$V2 <- sub("S", "E", sex$V2)
sex <- data.frame(t(sex[, c(2,5)]))
colnames(sex) <- sex[1, ]
sex <- sex[-1, ]
sex <- sex[, order(colnames(sex))]
sex <- sex[, which(names(sex) %in% names(pca))]
sex <- cbind(Covariates = "Sex", sex)

cov <- rbind(sex, pca[1:aha2, ])

write.table(cov, paste0("huvecImputeCov/", aha1, "_", "Sex", "_", aha2, "pc.txt"), sep = "\t", row.names = F, quote = F)

cat("Covariates file Ready  \n")
