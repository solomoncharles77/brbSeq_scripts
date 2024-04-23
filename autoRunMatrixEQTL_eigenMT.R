# Run MatrixEQTL and eigenMT
# Charles Solomon
# 15/07/2022
# Code credit -------------------------------------------------------------
# Most of the code on this script were written by Qin Qin Huang, and copied from 
# https://github.com/QinqinHuang/CAS_eQTL   and
# https://github.com/QinqinHuang/eQTL_simulations

# This script takes three arguments 
#  - the prefix name of the plink pca output
#  - number of genotype PC to use
#  - number of peer factors to use

# Assign command line arguments
aha <- commandArgs(trailingOnly = TRUE)
aha1 <- aha[1]
aha2 <- as.numeric(aha[2])
aha3 <- as.numeric(aha[3])

# load library
suppressPackageStartupMessages(library(MatrixEQTL))
suppressPackageStartupMessages(library(qqman))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doMC))

# files addresses
snpFile <- paste0("genoFiles/", aha1, "_qtlReady.txt")
geneFile <- paste0("exprFiles/", aha1, "_geneExpr_qtlReady.txt")
cvFile <- paste0("covFiles/", aha1, "_", aha2, "pc_", aha3,"peer.txt")
snpsLocation <- paste0("genoFiles/", aha1, ".Geno_snpsloc.txt")
genesLocation <- paste0("exprFiles/", aha1, "_genePos_qtlReady.txt")

# create output folders
resFilePrefix <- paste0(aha1, "_cisEQTL_", aha2, "pc_", aha3,  "pf")
resDir <- paste0(aha1, "_cisEQTL_", aha2, "pc_", aha3,  "pf_Results/")
if (dir.exists(resDir)) {
  unlink(resDir, recursive = TRUE)
  dir.create(resDir)
}else{
  dir.create(resDir)
}

eigenInDir <- paste0("eigenMT/eMT_input/input_QTL_", resFilePrefix, "/")
if (dir.exists(eigenInDir)) {
  unlink(eigenInDir, recursive = TRUE)
  dir.create(eigenInDir)
}else{
  dir.create(eigenInDir)
}

eigenOutDir <- paste0("eigenMT/eMT_output/output_QTL_", resFilePrefix, "/")
if (dir.exists(eigenOutDir)) {
  unlink(eigenOutDir, recursive = TRUE)
  dir.create(eigenOutDir)
}else{
  dir.create(eigenOutDir)
}

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(snpFile);

## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(geneFile);


## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$LoadFile(cvFile);


## import coordinate files
snpspos = read.table(snpsLocation, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(genesLocation, header = TRUE, stringsAsFactors = FALSE);

# meTest was run with NA filtered snp data with no covariates, no coordinates data.
me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name  = NULL,
  pvOutputThreshold  = 0,
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name.cis = NULL,
  pvOutputThreshold.cis = 0.09,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = 1000000,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);


save(me, file = paste0(resDir, resFilePrefix, ".RData"))

# QQplot
png(paste0(resDir, resFilePrefix, "_QQPlot.png"), width=6*300, height=6*300, res=300)
plot(me)
dev.off()

# Explore cis eQTL results -------------------------------------------------
colnames(snpspos) <- c("snps", "chr", "pos")

# subset cis data
cis <- me$cis$eqtls

# match join both datasets by snp
cisPos <- merge(snpspos, cis, by = "snps")

# png(paste0(resDir, resFilePrefix, "_Manhattan.png"), width=6*300, height=6*300, res=300)
# #manhattan plot
# manhattan(cisPos, chr = "chr", bp = "pos", p = "pvalue", snp = "snps",
#           col = c("gray10", "gray60"), chrlabs = NULL, logp = TRUE,
#           suggestiveline = FALSE, genomewideline = FALSE)
# dev.off()


# Export all raw cis eQTLs --------------------------------------------------------
fwrite(cis, paste0(resDir, resFilePrefix, ".csv"))
fwrite(cisPos, paste0(resDir, resFilePrefix, "_Position.csv"))

# # Export cis eQTLs by chromosome as required by eigenMT --------------------------------
# cis <- cis[, c("snps","gene","pvalue")]
# 
# ntests <- foreach(chr = 1:22, .combine = "+") %dopar% {
#   cat("  Saving MatrixeQTL results for Chr",chr,"...\n")
#   # SNPs on this chromosome
#   chrsnps <- snpspos[which(snpspos$chr == chr),]
#   # Tests involving these SNPs
#   eigenMTinput <- cis[which(cis$snps %in% chrsnps$snps),]
#   write.table(eigenMTinput, paste0(eigenInDir, "chr",chr,"_", resFilePrefix, "_4MT_pval.txt"),
#               quote = F, sep = "\t", row.names = F)
#   return(nrow(eigenMTinput))
# }
# 
# #system(paste0("python huvecImputeGenoEQTL_scripts/eigenMT.py --CHROM 19 --QTL ", eigenInDir, "chr",chr,"_", resFilePrefix, "_4MT_pval.txt --GEN eigenMT/inputGeno/",  aha1, "_SNP_chr19.txt --GENPOS eigenMT/inputGeno/",  aha1, "_snpsloc_chr19.txt --PHEPOS exprFiles/", aha1, "_genePos_qtlReady.txt --OUT ", eigenOutDir, "chr",chr,"_", resFilePrefix, "_eigenMT_FDR.txt"))
# 
# system(paste0("for chr in {1..22}
#               do
#               echo Processing CHR ${chr}
#               python huvecImputeGenoEQTL_scripts/eigenMT.py --CHROM ${chr} --QTL ", eigenInDir, "chr${chr}", "_", resFilePrefix, "_4MT_pval.txt --GEN eigenMT/inputGeno/",  aha1, "_SNP_chr${chr}.txt --GENPOS eigenMT/inputGeno/",  aha1, "_snpsloc_chr${chr}.txt --PHEPOS exprFiles/", aha1, "_genePos_qtlReady.txt --OUT ", eigenOutDir, "chr${chr}", "_", resFilePrefix, "_eigenMT_FDR.txt
#               done"))
