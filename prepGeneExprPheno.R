
# # This script takes one argument -  the prefix name of the plink processed output
# aha <- commandArgs(trailingOnly = TRUE)
# 
# aha2 <- paste0("genoFiles/", aha, ".Geno.txt")
# aha3 <- paste0("phenoFiles/", aha, ".txt")

# load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tximport))


# Import and process geno data --------------------------------------------
genoSamps <- read.table("genoFiles/brbSeqGeno_samples2.txt")


cat("\n")

cat("This is a view of the expression data \n")
cat("\n")

# Import and normalize count data --------------------------------------------
count <- read.csv("testSeq/readCount/rawReadCount.csv")
row.names(count) <- count$X
count$X <- NULL
count <- split.default(count, sub('\\d+', '', names(count)))


# Processing TNFA expression data -----------------------------------------
cat("\n")
cat("Processing TNFA expression data \n")
cat("\n")

tnfaCountdata <- count$X_TNFa
colnames(tnfaCountdata) <- sub("_TNFa", "", colnames(tnfaCountdata))
colnames(tnfaCountdata) <- sub("X", "S", colnames(tnfaCountdata))
setdiff(colnames(tnfaCountdata), genoSamps$V2)


cat("\n")
cat("Normalizing expression TNFA data \n")
cat("\n")

# normalize expr
sampleTable_TNFA <- data.frame(sampleID = colnames(tnfaCountdata), condition = factor("TNFA"))

dds_TNFA <- DESeqDataSetFromMatrix(countData = tnfaCountdata , colData = sampleTable_TNFA, ~1)
ddsF_TNFA <- dds_TNFA[ rowSums(counts(dds_TNFA)) > ncol(dds_TNFA)/2, ]
vst_TNFA=varianceStabilizingTransformation(ddsF_TNFA)
normExpr_TNFA <- as.data.frame(assay(vst_TNFA))

# Rename to match name in vcf file, order samples and Re-add geneID columns
normExpr_TNFA <- normExpr_TNFA[, order(names(normExpr_TNFA))]
normExpr_TNFA <- cbind(geneID = rownames(normExpr_TNFA), normExpr_TNFA)

cat("This is a view of the normalized TNFA expression data \n")
cat("\n")
normExpr_TNFA[1:10, 1:10]


# Prep bed for normalized genes -----------------------------------
# Add coordinate information 
exprCoord <- data.frame(fread("/scratch/vasccell/cs806/exprPhenoData/VSMC_Gene_Coords_withStrand.txt"))
normExprCoord_TNFA <- exprCoord[which(exprCoord$geneID %in% normExpr_TNFA$geneID), ]
normExprCoord_TNFA <- merge(normExprCoord_TNFA, normExpr_TNFA, by = "geneID")
normExprCoord_TNFA <- normExprCoord_TNFA[, c(2:4,1,1,5, 6:ncol(normExprCoord_TNFA))]
colnames(normExprCoord_TNFA)[1:6] <- c("#chr", "start", "end", "feature", "gene", "strand")
normExprCoord_TNFA <- normExprCoord_TNFA[order(normExprCoord_TNFA[,1], normExprCoord_TNFA[,2]), ]
normExprCoord_TNFA[1:10, 1:10]
cat("\n")
cat("Exporting QTL ready files \n")
cat("\n")

# Export files for QTLtools ---------------------------------------------
exprFileName_TNFA <- paste0("phenoFiles/normGeneExpr_TNFA_qtlTools_Ready.bed")
fwrite(normExprCoord_TNFA, file = exprFileName_TNFA, sep = "\t")

# system(paste0("module load samtools"))
# system(paste0("module load tabix"))
system(paste0("bgzip -f ", exprFileName_TNFA))
system(paste0("tabix ", exprFileName_TNFA, ".gz"))

###########################################################################################



# Processing UNTR expression data -----------------------------------------
cat("\n")
cat("Processing UNTR expression data \n")
cat("\n")

untrCountdata <- count$X_untr
colnames(untrCountdata) <- sub("_untr", "", colnames(untrCountdata))
colnames(untrCountdata) <- sub("X", "S", colnames(untrCountdata))
setdiff(colnames(untrCountdata), genoSamps$V2)


cat("\n")
cat("Normalizing expression UNTR data \n")
cat("\n")

# normalize expr
sampleTable_UNTR <- data.frame(sampleID = colnames(untrCountdata), condition = factor("UNTR"))

dds_UNTR <- DESeqDataSetFromMatrix(countData = untrCountdata , colData = sampleTable_UNTR, ~1)
ddsF_UNTR <- dds_UNTR[ rowSums(counts(dds_UNTR)) > ncol(dds_UNTR)/2, ]
vst_UNTR=varianceStabilizingTransformation(ddsF_UNTR)
normExpr_UNTR <- as.data.frame(assay(vst_UNTR))

# Rename to match name in vcf file, order samples and Re-add geneID columns
normExpr_UNTR <- normExpr_UNTR[, order(names(normExpr_UNTR))]
normExpr_UNTR <- cbind(geneID = rownames(normExpr_UNTR), normExpr_UNTR)

cat("This is a view of the normalized UNTR expression data \n")
cat("\n")
normExpr_UNTR[1:10, 1:10]


# Prep bed for normalized genes -----------------------------------
# Add coordinate information 
exprCoord <- data.frame(fread("/scratch/vasccell/cs806/exprPhenoData/VSMC_Gene_Coords_withStrand.txt"))
normExprCoord_UNTR <- exprCoord[which(exprCoord$geneID %in% normExpr_UNTR$geneID), ]
normExprCoord_UNTR <- merge(normExprCoord_UNTR, normExpr_UNTR, by = "geneID")
normExprCoord_UNTR <- normExprCoord_UNTR[, c(2:4,1,1,5, 6:ncol(normExprCoord_UNTR))]
colnames(normExprCoord_UNTR)[1:6] <- c("#chr", "start", "end", "feature", "gene", "strand")
normExprCoord_UNTR <- normExprCoord_UNTR[order(normExprCoord_UNTR[,1], normExprCoord_UNTR[,2]), ]
normExprCoord_UNTR[1:10, 1:10]
cat("\n")
cat("Exporting QTL ready files \n")
cat("\n")

# Export files for QTLtools ---------------------------------------------
exprFileName_UNTR <- paste0("phenoFiles/normGeneExpr_UNTR_qtlTools_Ready.bed")
fwrite(normExprCoord_UNTR, file = exprFileName_UNTR, sep = "\t")

# system(paste0("module load samtools"))
# system(paste0("module load tabix"))
system(paste0("bgzip -f ", exprFileName_UNTR))
system(paste0("tabix ", exprFileName_UNTR, ".gz"))

###########################################################################################
