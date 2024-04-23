
# This script takes one argument -  the prefix name of the plink processed output
aha <- commandArgs(trailingOnly = TRUE)

aha2 <- paste0("genoFiles/", aha, ".Geno.txt")
aha3 <- paste0("exprFiles/", aha, ".txt")

# load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))


# Import and process geno data --------------------------------------------
geno <- data.frame(fread(aha2))
invisible(gc())
colnames(geno) <- gsub(".*_S", "S", colnames(geno))
rownames(geno) <- geno$SNP
geno$SNP <- NULL

cat("\n")

cat("This is a view of the expression data \n")
cat("\n")
# Import and process expr data --------------------------------------------
expr <- data.frame(fread(aha3, header = TRUE))
# Harmonize names and samples in geno and expr ----------------------------
expr[1:10, 1:10]

# Extract and normalize expr samples present in geno data
expr1 <- expr[, which(colnames(expr) %in% colnames(geno))]
rownames(expr1) <- expr$GeneID


cat("\n")

cat("This is a view of the genotype data \n")
cat("\n")

geno <- geno[, order(colnames(geno))] 
geno[1:10, 1:10]

cat("\n")
cat("Normalizing expression data \n")
cat("\n")

# normalize expr
sampleTable <- data.frame(condition = factor(rep(c("HUVEC"), ncol(expr1))))
rownames(sampleTable) <- colnames(expr1)

dds <- DESeqDataSetFromMatrix(round(expr1), sampleTable, ~1)
ddsF <- dds[ rowSums(counts(dds)) > ncol(dds), ]
vst=varianceStabilizingTransformation(ddsF)
normExpr <- as.data.frame(assay(vst))


cat("This is a view of the normalized expression data \n")
cat("\n")

normExpr[1:10, 1:10]

# Check if samples match in geno and expr and Re-add gene/SNP-ID columns
normExpr <- normExpr[, order(names(normExpr))]
geno <- geno[, order(names(geno))]
all(colnames(geno) == colnames(normExpr))

normExpr <- cbind(geneID = rownames(normExpr), normExpr)
geno <- cbind(SNP = rownames(geno), geno)

# Extract location for normalized genes -----------------------------------
exprCoord <- data.frame(fread("/scratch/vasccell/cs806/exprPhenoData/VSMC_Gene_Coords.txt"))
normExprCoord <- exprCoord[which(exprCoord$geneID %in% rownames(normExpr)), ]

cat("\n")
cat("Exporting eQTL ready files \n")
cat("\n")

# Export files for matrixEQTL ---------------------------------------------
genoFileName <- paste0("genoFiles/", aha, "_qtlReady.txt")
fwrite(geno, file = genoFileName, sep = "\t")

exprFileName <- paste0("exprFiles/", aha, "_geneExpr_qtlReady.txt")
fwrite(normExpr, file = exprFileName, sep = "\t")

coordFileName <- paste0("exprFiles/", aha, "_genePos_qtlReady.txt")
fwrite(normExprCoord, file = coordFileName, sep = "\t", quote = F)
