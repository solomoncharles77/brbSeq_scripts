
library(data.table)
library(tidyverse)
library(ggvenn)
library(openxlsx)

# qtlHeader <- read.delim("brbSeq_scripts/QTLtools_nominal_header.txt", header = F)
# colIDs <- qtlHeader$V2
# colIDs <- sub(" \\| ve_by_pc1 \\| n_phe_in_grp", "", colIDs)
# colIDs <- sub("phe_id \\| ", "", colIDs)
# 
# tnfa <- data.frame(fread("geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL/geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL_nominal_ALL_pval0.05.txt.gz"))
# colnames(tnfa) <- paste0("tnfa_", colIDs)
# 
# untr <- data.frame(fread("geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL/geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL_nominal_ALL_pval0.05.txt.gz"))
# colnames(untr) <- paste0("untr_", colIDs)
# gc()
# 
# tnfaUntr <- merge(tnfa, untr, by.x = c("tnfa_phe_id", "tnfa_var_id"), by.y = c("untr_phe_id", "untr_var_id"))
# 
# head(tnfa)
# 
# 
# # Detecting reQTLs by eQTL β-comparison --------------------------------------------------------------
# 
# tnfaUntr$zScore <- (tnfaUntr$untr_slope - tnfaUntr$tnfa_slope) / sqrt(tnfaUntr$untr_slope_se^2 + tnfaUntr$tnfa_slope_se^2)
# #tnfaUntr$zScorePval <- pnorm(tnfaUntr$zScore, mean = 0, sd = 1, lower.tail = TRUE)
# tnfaUntr$zScorePval <- 2 * pnorm(-abs(tnfaUntr$zScore))
# tnfaUntr2 <- tnfaUntr[tnfaUntr$zScorePval < 0.05, ]


###############################################################################

tnfa <- data.frame(fread("geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL/geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL_permute_ALL_rsID.txt.gz"))
colnames(tnfa) <- paste0("tnfa_", colnames(tnfa))
tnfa <- na.omit(tnfa)
write.csv(tnfa, "resFiles/geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL_permute_ALL_annot.csv", row.names = F)
tnfa <- tnfa[tnfa$tnfa_adj_emp_pval < 0.05, ]
write.csv(tnfa, "resFiles/geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL_permute_Sig_annot.csv", row.names = F)
write.table(tnfa$tnfa_grp_id, "resFiles/geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL_candGenes.txt", row.names = F, col.names = F, quote = F)
write.table(tnfa[, c(1,30,2,24)], "resFiles/geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL_topAssocs.txt", sep = "\t", row.names = F, col.names = F, quote = F)


untr <- data.frame(fread("geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL/geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL_permute_ALL_rsID.txt.gz"))
colnames(untr) <- paste0("untr_", colnames(untr))
untr <- na.omit(untr)
write.csv(untr, "resFiles/geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL_permute_ALL_annot.csv", row.names = F)
untr <- untr[untr$untr_adj_emp_pval < 0.05, ]
write.csv(untr, "resFiles/geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL_permute_Sig_annot.csv", row.names = F)
write.table(untr$untr_grp_id, "resFiles/geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL_candGenes.txt", row.names = F, col.names = F, quote = F)
write.table(untr[, c(1,30,2,24)], "resFiles/geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL_topAssocs.txt", sep = "\t", row.names = F, col.names = F, quote = F)


ggplot(data = tnfa, aes(y = tnfa_slope, x = tnfa_maf)) +
  geom_point() +
  labs(title = "TNFa eQTL Beta vs MAF Plot",
       x = "Minor Allele Frequency",
       y = "Effect Size") +
  theme_minimal() -> tnfaEqtlMafBetaPlot
tnfaEqtlMafBetaPlot
ggsave("tnfaEqtlMafBetaPlot.png", plot = tnfaEqtlMafBetaPlot,
       path = "resPlots", height = 4.73, width = 7.45, units = "in", bg = "white")

ggplot(data = untr, aes(y = untr_slope, x = untr_maf)) +
  geom_point() +
  labs(title = "Naive eQTL Beta vs MAF Plot",
       x = "Minor Allele Frequency",
       y = "Effect Size") +
  theme_minimal() -> naiveEqtlMafBetaPlot
naiveEqtlMafBetaPlot
ggsave("naiveEqtlMafBetaPlot.png", plot = naiveEqtlMafBetaPlot,
       path = "resPlots", height = 4.73, width = 7.45, units = "in", bg = "white")


# Get eGene overlap -------------------------------------------------------

tnfa2 <- tnfa[!duplicated(tnfa$tnfa_phe_id), ]
untr2 <- untr[!duplicated(untr$untr_phe_id), ]
ovList <- list(untr2$untr_phe_id, tnfa2$tnfa_phe_id)
names(ovList) <- c("NAIVE", "TNFA")
tnfaUntrOverlap <- gplots::venn(ovList)

ggplot() + 
  geom_venn(ovList, size=2, textsize=2.5) +
  scale_fill_manual(values=c("red","blue")) +
  theme_void() +
  theme(legend.position = "none")  -> ovVenn

ggsave("naiveTNFA_eGene_overlap.pdf", plot = ovVenn,
       path = "resPlots/Fig2", height = 1.5, width = 1.5, units = "in", bg = "white")

write.xlsx(tnfaUntrOverlap, "resFiles/naive_tnfa_eGene_overlap.xlsx")

sharedTNFA <- tnfa2[tnfa2$ensembl_gene_id %in% attr(tnfaUntrOverlap,"intersections")$`NAIVE:TNFA`, ]
sharedUNTR <- untr2[untr2$ensembl_gene_id %in% attr(tnfaUntrOverlap,"intersections")$`NAIVE:TNFA`, ]

sharedTnfaUntr <- merge(sharedTNFA, sharedUNTR, by.x = c("tnfa_phe_id", "tnfa_var_id"), by.y = c("untr_phe_id", "untr_var_id"))


##################################################################################
##################################################################################
# sQTL

tnfa <- data.frame(fread("leafcutterSplicing_TNFA_cisEQTL/leafcutterSplicing_TNFA_cisEQTL_permute_ALL.txt.gz"))
colnames(tnfa) <- paste0("tnfa_", colIDs)
tnfa <- merge(tnfa[, c(1,6,9,10,18,20:23)], maf[, c(1, 4:6)], by.x = "tnfa_var_id", by.y = "var_id")
tnfa <- merge(annot[, -3], tnfa, by.x = "ensembl_gene_id", by.y = "tnfa_grp_id")
tnfa <- na.omit(tnfa)
write.csv(tnfa, "resFiles/leafcutterSplicing_TNFA_cisEQTL_permute_ALL_annot.csv", row.names = F)
tnfa <- tnfa[tnfa$tnfa_adj_emp_pval < 0.05, ]
write.table(tnfa$ensembl_gene_id, "resFiles/leafcutterSplicing_TNFA_cisSQTL_candGenes.txt", row.names = F, col.names = F, quote = F)


untr <- data.frame(fread("leafcutterSplicing_UNTR_cisEQTL/leafcutterSplicing_UNTR_cisEQTL_permute_ALL.txt.gz"))
colnames(untr) <- paste0("untr_", colIDs)
untr <- merge(untr[, c(1,6,9,10,18,20:23)], maf[, c(1, 4:6)], by.x = "untr_var_id", by.y = "var_id")
untr <- merge(annot[, -3], untr, by.x = "ensembl_gene_id", by.y = "untr_grp_id")
untr <- na.omit(untr)
write.csv(untr, "resFiles/leafcutterSplicing_NAIVE_cisEQTL_permute_ALL_annot.csv", row.names = F)
untr <- untr[untr$untr_adj_emp_pval < 0.05, ]
write.table(untr$ensembl_gene_id, "resFiles/leafcutterSplicing_UNTR_cisSQTL_candGenes.txt", row.names = F, col.names = F, quote = F)


ggplot(data = tnfa, aes(y = tnfa_slope, x = maf)) +
  geom_point() +
  labs(title = "TNFa sQTL Beta vs MAF Plot",
       x = "Minor Allele Frequency",
       y = "Effect Size") +
  theme_minimal() -> tnfaSqtlMafBetaPlot
tnfaSqtlMafBetaPlot
ggsave("tnfaSqtlMafBetaPlot.png", plot = tnfaSqtlMafBetaPlot,
       path = "resPlots", height = 4.73, width = 7.45, units = "in", bg = "white")

ggplot(data = untr, aes(y = untr_slope, x = maf)) +
  geom_point() +
  labs(title = "Naive sQTL Beta vs MAF Plot",
       x = "Minor Allele Frequency",
       y = "Effect Size") +
  theme_minimal() -> naiveSqtlMafBetaPlot
naiveSqtlMafBetaPlot
ggsave("naiveSqtlMafBetaPlot.png", plot = naiveSqtlMafBetaPlot,
       path = "resPlots", height = 4.73, width = 7.45, units = "in", bg = "white")

# Get eGene overlap -------------------------------------------------------

tnfa2 <- tnfa[!duplicated(tnfa$tnfa_phe_id), ]
untr2 <- untr[!duplicated(untr$untr_phe_id), ]
ovList <- list(untr2$untr_phe_id, tnfa2$tnfa_phe_id)
names(ovList) <- c("NAIVE", "TNFA")
tnfaUntrOverlap <- gplots::venn(ovList)

tnfa2 <- tnfa[!duplicated(tnfa$ensembl_gene_id), ]
untr2 <- untr[!duplicated(untr$ensembl_gene_id), ]
ovList <- list(untr2$ensembl_gene_id, tnfa2$ensembl_gene_id)
names(ovList) <- c("NAIVE", "TNFA")
tnfaUntrOverlap <- gplots::venn(ovList)


ggplot() + 
  geom_venn(ovList, size=4, textsize=5) +
  scale_fill_manual(values=c("red","blue")) +
  theme_void() +
  theme(legend.position = "none") -> ovVenn

ggsave("naiveTNFA_sGene_overlap.png", plot = ovVenn,
       path = "resPlots", height = 4.73, width = 7.45, units = "in", bg = "white")

write.xlsx(tnfaUntrOverlap, "resFiles/naive_tnfa_sGene_overlap.xlsx")


