suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(tidyverse))

# Import data -------------------------------------------------------------
# expr1 <- data.frame(fread("phenoFiles/normGeneExpr_UNTR_qtlTools_Ready.bed.gz"))
# expr2 <- data.frame(fread("phenoFiles/normGeneExpr_TNFA_qtlTools_Ready.bed.gz"))

expr1 <- data.frame(fread("phenoFiles/geneExpr_TMM_IVTnorm_UNTR_qtlTools_Ready.bed.gz"))
expr2 <- data.frame(fread("phenoFiles/geneExpr_TMM_IVTnorm_TNFA_qtlTools_Ready.bed.gz"))

expr1 <- expr1[, c(5, 7:ncol(expr1))]
expr2 <- expr2[, c(5, 7:ncol(expr2))]
expr3 <- merge(expr1, expr2, by = "gene")
colnames(expr3) <- gsub("\\..*", "", colnames(expr3))


# cov1 <- data.frame(fread("covFiles/normGeneExpr_UNTR_geno3pc_pheno40pc.txt"))
# cov2 <- data.frame(fread("covFiles/normGeneExpr_TNFA_geno3pc_pheno40pc.txt"))

cov1 <- data.frame(fread("covFiles/geneExpr_TMM_IVTnorm_UNTR_geno3pc_pheno40pc.txt"))
cov2 <- data.frame(fread("covFiles/geneExpr_TMM_IVTnorm_TNFA_geno3pc_pheno40pc.txt"))


cov3 <- cbind(cov1, cov2[, -c(1)])
cov3$SampleID <- sub("covFiles/normGeneExpr_UNTR_qtlTools_Ready_1_1_svd_", "pheno", cov3$SampleID)
rownames(cov3) <- cov3$SampleID
cov3 <- data.frame(t(cov3[,-1]))


geno <- data.frame(fread("genoFiles/brbSeqGeno_QCed_Geno.txt"))
colnames(geno) <- gsub(".*_S", "S", colnames(geno))
geno2 <- cbind(geno, geno[, -c(1)])


compareQTLs <-function(cands, expr3, geno2, cov3) {
  # Register parallel backend to use multiple cores
  registerDoParallel(cores = 14)  
  
  # Execute the foreach loop
  results <- foreach(rr = 1:nrow(cands), .combine = rbind) %dopar% {
    # Extract relevant data from the inputs
    candGene <- cands$geneID[rr]
    geneName <- cands$geneName[rr]
    candSNP <- cands$var_id[rr]
    
    # Filter expression and genotype data
    candExpr <- expr3[expr3$gene %in% candGene, ]
    candGeno <- geno2[geno2$SNP %in% candSNP, ]
    
    # Merge and transpose data
    df <- t(rbind(candGeno[-1], candExpr[-1]))
    if (ncol(df) == 1) {
      cat(candGene, " is absent from merged expression data \n")
      return(NULL)  # Skip to the next iteration
    }
    colnames(df) <- c("snp", "gene_expr")
    df <- data.frame(df)
    df2 <- merge(df, cov3, by=0)
    df2$condition <- as.factor(ifelse(grepl(".1", df2$Row.names), 1, 0))
    df2$Donor <- as.factor(gsub("\\.1", "", df2$Row.names))
    df2$Row.names <- NULL
    df2 <- df2[complete.cases(df2), ]
    
    # Create and update model formulas
    covIDs <- colnames(cov3)
    mf <- reformulate(covIDs, response = "gene_expr")
    mf1 <- update(mf, . ~ . + snp + condition + (1 | Donor))
    mf2 <- update(mf, . ~ . + snp * condition + (1 | Donor))
    
    # Fit models and compare
    h0 <- lmer(mf1, data = df2, REML=F)
    h1 <- lmer(mf2, data = df2, REML=F)
    anovaRes <- anova(h0, h1)
    anovaRes <- broom::tidy(anovaRes)
    anovaRes$geneID <- candGene
    anovaRes$geneName <- geneName
    anovaRes$var_id <- candSNP
    
    return(anovaRes)
  }
  
  # Stop parallel backend
  stopImplicitCluster()
  
  return(results)
}

untrCands <- read.csv("resFiles/geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL_permute_Sig_annot.csv")
untrCands <- setnames(untrCands, c("untr_grp_id", "untr_external_gene_name", "untr_var_id"),
                      c("geneID", "geneName", "var_id"))
untrSpec <- compareQTLs(untrCands, expr3, geno2, cov3)
untrSpec <- untrSpec[untrSpec$term != "h0", ]
untrSpecSig <- untrSpec[untrSpec$p.value < 0.05, ]
untrSpecSig <- untrCands[untrCands$geneID %in% untrSpecSig$geneID, ]
write.csv(untrSpecSig, "resFiles/naive_condition_eQTL.csv", row.names = F)
write.table(untrSpecSig[, c(1,30,2,24)], "resFiles/ncsCands.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(untrSpecSig[, c(1)], "resFiles/ncsCands_geneList.txt", sep = "\t", quote = F, row.names = F, col.names = F)


tnfaCands <- read.csv("resFiles/geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL_permute_Sig_annot.csv")
tnfaCands <- setnames(tnfaCands, c("tnfa_grp_id", "tnfa_external_gene_name", "tnfa_var_id"),
                      c("geneID", "geneName", "var_id"))
tnfaSpec <- compareQTLs(tnfaCands, expr3, geno2, cov3)
tnfaSpec <- tnfaSpec[tnfaSpec$term != "h0", ]
tnfaSpecSig <- tnfaSpec[tnfaSpec$p.value < 0.05, ]
tnfaSpecSig <- tnfaCands[tnfaCands$geneID %in% tnfaSpecSig$geneID, ]
write.csv(tnfaSpecSig, "resFiles/tnfa_response_eQTL.csv", row.names = F)
write.table(tnfaSpecSig[, c(1,30,2,24)], "resFiles/trsCands.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(tnfaSpecSig[, c(1)], "resFiles/trsCands_geneList.txt", sep = "\t", quote = F, row.names = F, col.names = F)


ovList <- list(untrSpecSig$geneID, tnfaSpecSig$geneID)
names(ovList) <- c("NAIVE", "TNFA")
tnfaUntrOverlap <- gplots::venn(ovList)

ggplot() + 
  geom_venn(ovList, size=4, textsize=5) +
  scale_fill_manual(values=c("red","blue")) +
  theme_void() +
  theme(legend.position = "none") -> ovVenn

ggsave("naiveTNFA_specSig_overlap.png", plot = ovVenn,
       path = "resPlots", height = 4.73, width = 7.45, units = "in", bg = "white")


