# Multiple testing correction of eQTL result
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

# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(org.Hs.eg.db))

# files addresses
resFilePrefix <- paste0(aha1, "_cisEQTL_", aha2, "pc_", aha3,  "pf")
resDir <- paste0(aha1, "_cisEQTL_", aha2, "pc_", aha3,  "pf_Results/")
eigenOutDir <- paste0("eigenMT/eMT_output/output_QTL_", resFilePrefix, "/")

# create output folders
csvDir <- paste0(resDir, resFilePrefix, "_SigAssocCSV/")
if (dir.exists(csvDir)) {
  unlink(csvDir, recursive = TRUE)
  dir.create(csvDir)
}else{
  dir.create(csvDir)
}

txtDir <- paste0(resDir, resFilePrefix, "_SigAssocTXT/")
if (dir.exists(txtDir)) {
  unlink(txtDir, recursive = TRUE)
  dir.create(txtDir)
}else{
  dir.create(txtDir)
}

topDir <- paste0(resDir, resFilePrefix, "_SigAssocTOP/")
if (dir.exists(topDir)) {
  unlink(topDir, recursive = TRUE)
  dir.create(topDir)
}else{
  dir.create(topDir)
}


# Functions for multiple testing correction -----------------------------------------
getThreshold <- function(d, opt_fdr = 0.05) {
  set0 <- d[which(d[,2] <= opt_fdr),]
  set1 <- d[which(d[,2] > opt_fdr),]
  threshold <- (sort(set1[,1])[1] - sort(-1.0 * set0[,1])[1]) / 2
  return(threshold)
}


getAssociation.withThresholds <- function(d1, d2, opt_fdr = 0.05) {
  d1 <- d1[which(d1$pvalue < opt_fdr),]
  d1 <- merge(d1, d2, by = "gene")
  d1 <- d1[which(d1$pvalue <= d1[,4]),]
  d1 <- d1[order(d1$pvalue), -4]
  return(d1)
}


mul_corr_eigenMT <- function(alltests, eMT_output, fdr_desired = 0.05, step2_met = "BH") {
  cat("  ** Step 1: eigenMT with window size 200, Step 2: ", step2_met, "-FDR ** \n", sep = "")
  # Keep the best association for each gene
  eMT_output <- merge(eMT_output, alltests[which(!duplicated(alltests$gene)), c("gene","pvalue")], by = "gene")
  # Calculate eMT adjusted pvalues for each top association
  eMT_output$eMT <- eMT_output$pvalue * eMT_output$TESTS
  eMT_output$eMT[which(eMT_output$eMT > 1)] <- 1
  
  # Correct for multiple genes
  if(step2_met == "BH") {
    eMT_output$gene.fdr <- p.adjust(eMT_output$eMT, "fdr")
  } else if(step2_met == "ST") {
    eMT_output$gene.fdr <- qvalue(eMT_output$eMT)$qvalues
    cat("    calculating qvalues in Step 2, pi0 =", qvalue(eMT_output$eMT)$pi0, "\n")
  } else if(step2_met == "BY") {
    eMT_output$gene.fdr <- p.adjust(eMT_output$eMT, "BY")
  } else if(step2_met == "Bonferroni") {
    eMT_output$gene.fdr <- p.adjust(eMT_output$eMT, "bonferroni")
  }
  
  # Get the empirical p-value threshold
  threshold <- getThreshold(eMT_output[,c("eMT","gene.fdr")], opt_fdr = fdr_desired)
  cat("    The empirical p-value threshold =", threshold, "\n")
  # Calculate the nominal pvalue threshold for each gene
  nthre.eMT <- data.frame(gene = eMT_output$gene, eMT = threshold / eMT_output$TESTS)
  names(nthre.eMT)[2] <- paste0("eigenMT.", step2_met)
  # Get the significant pairs based on the nominal thresholds
  dd <- getAssociation.withThresholds(d1 = alltests[, c("snps","gene","pvalue")], d2 = nthre.eMT)
  names(dd)[3] <- paste0("eigenMT.", step2_met)
  
  return(list(dd,nthre.eMT))
}

################################################################################################


# Read MatrixeQTL cis eQTL output
cis <- fread(paste0(resDir, resFilePrefix, ".csv"))
# Sort by absulute statistics
cis <- cis[order(abs(cis$statistic), decreasing = T),]


# Do MT at once in a loop -------------------------------------------------
ntests <- foreach(chr = 1:22, .combine = "+") %dopar% {
  cat("  reading eigenMT output for Chr",chr,"...\n")
  dd <- read.table(paste0(eigenOutDir, "chr",chr,"_", resFilePrefix, "_eigenMT_FDR.txt"), header = T)
  dd1 <- dd[which(complete.cases(dd$TESTS)), c("gene","TESTS")]
  cat("  performing eigenMT-BH multiple testing for Chr",chr,"...\n")
  #---------- eigenMT-BH ---------
  eMTBH_returnlist <- mul_corr_eigenMT(alltests = cis, eMT_output = dd1, step2_met = "BH")
  sigeAsso_eigenMTBH <- eMTBH_returnlist[[1]]
  # Nominal thresholds
  nominalthresholds <- eMTBH_returnlist[[2]]
  
  
  # Save results
  cat("  saving eigenMT-BH multiple testing for Chr",chr,"...\n")
  #write.table(sigeAsso_eigenMTBH, "sigeAsso_eigenMTBH.txt", quote = F, sep = "\t", row.names = F)
  #write.table(nominalthresholds, "Nominal_thresholds.txt", quote = F, sep = "\t", row.names = F)
  
  cis_sub <- cis[which(cis$gene %in% sigeAsso_eigenMTBH$gene & cis$snps %in% sigeAsso_eigenMTBH$snps), c("snps","gene","statistic","pvalue","beta")]
  dd <- merge(sigeAsso_eigenMTBH[, 1:2], cis_sub, by = c("gene","snps"))
  dd <- dd[order(abs(dd$statistic), decreasing = T), ]
  dd <- dd[, c( "gene", "snps", "statistic", "pvalue", "beta")]
  appel <- as.character(dd$gene)
  dd$gene_Name <- mapIds(org.Hs.eg.db,
                         keys=appel,
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")
  dd$snps <- sub("chr", "", dd$snps)
  snpList <- data.frame(snps = dd$snps)
  snpListFile <- paste0("getRS.txt")
  write.table(snpList, file = snpListFile, quote = F, row.names = F, col.names = F)
  # get rsID
  snpRS <- system(paste0("/home/c/cs806/tsv-utils-v2.1.2_linux-x86_64_ldc2/bin/tsv-join -f ", snpListFile, " -k 1 -d 2 /scratch/vasccell/cs806/colocalization/dbSNP/dbSNP155.hg38.rsID_CoordID.txt"), intern = T)
  snpRS <- data.frame(fread(text = snpRS, header = F))
  colnames(snpRS) <- c("rsID", "snps")
  dd <- merge(dd, snpRS, by = "snps", all.x = T)
  dd <- dd[, c( "gene", "gene_Name", "snps", "rsID", "statistic", "pvalue", "beta")]
  dd2 <- dd[which(!duplicated(dd$gene)), ]
  write.table(dd, paste0(txtDir, resFilePrefix,"_Chr", chr, "_sigeAsso_eigenMT-BH.txt"), quote = F, sep = "\t", row.names = F)  
  write.csv(dd, paste0(csvDir, resFilePrefix,"_Chr", chr, "_sigeAsso_eigenMT-BH.csv"), row.names = F)
  write.csv(dd2, paste0(topDir, resFilePrefix,"_Chr", chr, "_sigeAsso_eigenMT-BH.csv"), row.names = F)
  gc()
}


csvFiles <- list.files(csvDir)
allCSV <- lapply(csvFiles, function(x){
  cc <- read.csv(paste0(csvDir, x))
})
allCSV <- do.call(rbind, allCSV)
write.csv(allCSV, paste0(csvDir, resFilePrefix, "_All_sigeAsso_eigenMT-BH.csv"), row.names = F)


txtFiles <- list.files(txtDir)
alltxt <- lapply(txtFiles, function(x){
  cc <- read.table(paste0(txtDir, x), header = T)
})
alltxt <- do.call(rbind, alltxt)
write.table(alltxt, paste0(txtDir, resFilePrefix, "_All_sigeAsso_eigenMT-BH.txt"), quote = F, sep = "\t", row.names = F)

