
library(data.table)
library(tidyverse)

# qtlHeader <- read.delim("brbSeq_scripts/QTLtools_nominal_header.txt", header = F)
# nomIDs <- qtlHeader$V2
# nomIDs <- sub(" \\| ve_by_pc1 \\| n_phe_in_grp", "", nomIDs)
# nomIDs <- sub("phe_id \\| ", "", nomIDs)

qtlHeader <- read.delim("brbSeq_scripts/QTLtools_permute_header.txt", header = F)
pemIDs <- qtlHeader$V2
pemIDs <- sub(" \\| ve_by_pc1 \\| n_phe_in_grp", "", pemIDs)
pemIDs <- sub("phe_id \\| ", "", pemIDs)


# Compare TNFA nominal with UNTR perm -------------------------------------
tnfaNom <- data.frame(fread("geneExprTNFA_brbSeq_cisEQTL/geneExprTNFA_brbSeq_cisEQTL_nominal_ALL_pval0.05_rsID.txt.gz",
                         select=c(1,2,14,16,17,19,24,25)))
colnames(tnfaNom) <- paste0("tnfa_", colnames(tnfaNom))
gc()
          
# Import UNTR permute 
untrPem <- data.frame(fread("geneExprUNTR_brbSeq_cisEQTL/geneExprUNTR_brbSeq_cisEQTL_permute_ALL.txt.gz",
                            select = c(1,6,10,18,20:22)))
untrPem <- untrPem[complete.cases(untrPem), ]
colnames(untrPem) <- paste0("untr_", pemIDs[c(1,6,10,18,20:22)])
untrPem <- untrPem[untrPem$untr_adj_emp_pval < 0.05, ]

tnfaUntr <- merge(tnfaNom, untrPem, by.x = c("tnfa_grp_id", "tnfa_var_id"), by.y = c("untr_phe_id", "untr_var_id"))

tnfaUntr$zScore <- (tnfaUntr$untr_slope - tnfaUntr$tnfa_slope) / sqrt(tnfaUntr$untr_slope_se^2 + tnfaUntr$tnfa_slope_se^2)
#tnfaUntr$zScorePval <- pnorm(tnfaUntr$zScore, mean = 0, sd = 1, lower.tail = TRUE)
tnfaUntr$zScorePval <- 2 * pnorm(-abs(tnfaUntr$zScore))


###########################################################################################

untrNom <- data.frame(fread("geneExprUNTR_brbSeq_cisEQTL/geneExprUNTR_brbSeq_cisEQTL_nominal_ALL_pval0.05_rsID.txt.gz",
                            select=c(1,2,14,16,17,19,24,25)))
colnames(untrNom) <- paste0("untr_", colnames(untrNom))
gc()

tnfaPem <- data.frame(fread("geneExprTNFA_brbSeq_cisEQTL/geneExprTNFA_brbSeq_cisEQTL_permute_ALL.txt.gz",
                            select = c(1,6,10,18,20:22)))
tnfaPem <- tnfaPem[complete.cases(tnfaPem), ]
colnames(tnfaPem) <- paste0("tnfa_", pemIDs[c(1,6,10,18,20:22)])
tnfaPem <- tnfaPem[tnfaPem$tnfa_adj_emp_pval < 0.05, ]

untrTnfa <- merge(untrNom, tnfaPem, by.x = c("untr_grp_id", "untr_var_id"), by.y = c("tnfa_phe_id", "tnfa_var_id"))

untrTnfa$zScore <- (untrTnfa$untr_slope - untrTnfa$tnfa_slope) / sqrt(untrTnfa$untr_slope_se^2 + untrTnfa$tnfa_slope_se^2)
#untrTnfa$zScorePval <- pnorm(untrTnfa$zScore, mean = 0, sd = 1, lower.tail = TRUE)
untrTnfa$zScorePval <- 2 * pnorm(-abs(untrTnfa$zScore))


# Get shared --------------------------------------------------------------
shared <- intersect( tnfaUntr$tnfa_grp_id, untrTnfa$untr_grp_id )

untrTnfa_Shared <- untrTnfa[untrTnfa$untr_grp_id %in% shared, ]
untrTnfa_Shared_similar <- untrTnfa_Shared[untrTnfa_Shared$zScorePval > 0.05, ]
untrTnfa_Shared_contra <- untrTnfa_Shared[untrTnfa_Shared$zScorePval < 0.05, ]

tnfaUntr_Shared <- tnfaUntr[tnfaUntr$tnfa_grp_id %in% shared, ]
tnfaUntr_Shared_similar <- tnfaUntr_Shared[tnfaUntr_Shared$zScorePval > 0.05, ]
tnfaUntr_Shared_contra <- tnfaUntr_Shared[tnfaUntr_Shared$zScorePval < 0.05, ]

similar <- intersect(untrTnfa_Shared_similar$untr_grp_id, tnfaUntr_Shared_similar$tnfa_grp_id)
sharedSimilar <- untrTnfa_Shared[untrTnfa_Shared$untr_grp_id %in% similar, c(1,8,2,6) ]
write.csv(sharedSimilar, "resFiles/sharedSimilar_assocs.csv", row.names = F)
write.table(sharedSimilar, "resFiles/sharedSimilar_assocs.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(sharedSimilar[, c(1)], "resFiles/sharedSimilar_assocs_eGenes.txt", row.names = F, col.names = F, quote = F, sep = "\t")


sharedSimilar2 <- tnfaUntr_Shared[tnfaUntr_Shared$tnfa_grp_id %in% similar, c(1,8,2,6) ]
write.csv(sharedSimilar2, "resFiles/sharedSimilar_assocs2.csv", row.names = F)
write.table(sharedSimilar2, "resFiles/sharedSimilar_assocs2.txt", row.names = F, col.names = F, quote = F, sep = "\t")


contra <- intersect(untrTnfa_Shared_contra$untr_grp_id, tnfaUntr_Shared_contra$tnfa_grp_id)
sharedContra <- untrTnfa_Shared[untrTnfa_Shared$untr_grp_id %in% contra, c(1,8,2,6) ]
write.csv(sharedContra, "resFiles/sharedContra_assocs.csv", row.names = F)
write.table(sharedContra, "resFiles/sharedContra_assocs.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(sharedContra[, c(1)], "resFiles/sharedContra_assocs_eGenes.txt", row.names = F, col.names = F, quote = F, sep = "\t")


sharedContra2 <- tnfaUntr_Shared[tnfaUntr_Shared$tnfa_grp_id %in% contra, c(1,8,2,6) ]
write.csv(sharedContra2, "resFiles/sharedContra_assocs2.csv", row.names = F)
write.table(sharedContra2, "resFiles/sharedContra_assocs2.txt", row.names = F, col.names = F, quote = F, sep = "\t")


# Get naive Significant
naiveSignificant <- tnfaUntr[!tnfaUntr$tnfa_grp_id %in% shared, ]
naiveSignificant <- tnfaUntr[tnfaUntr$zScorePval < 0.05, ]
naiveSignificant <- naiveSignificant[, c(1,8,2,6,11,12,10,13) ]
colnames(naiveSignificant) <- sub("tnfa_", "", colnames(naiveSignificant))
colnames(naiveSignificant) <- sub("untr_", "", colnames(naiveSignificant))
colnames(naiveSignificant)[1:2] <- c("geneID", "geneName") 
write.csv(naiveSignificant, "resFiles/naive-Significant_eGenes.csv", row.names = F)
write.table(naiveSignificant[, c(1:4)], "resFiles/naive_Significant_eGenes.txt", row.names = F, col.names = F, quote = F, sep = "\t")


# Get tnfa Significant
tnfaSignificant <- untrTnfa[!untrTnfa$untr_grp_id %in% shared, ]
tnfaSignificant <- untrTnfa[untrTnfa$zScorePval < 0.05, ]
tnfaSignificant <- tnfaSignificant[, c(1,8,2,6,11,12,10,13) ]
colnames(tnfaSignificant) <- sub("tnfa_", "", colnames(tnfaSignificant))
colnames(tnfaSignificant) <- sub("untr_", "", colnames(tnfaSignificant))
colnames(tnfaSignificant)[1:2] <- c("geneID", "geneName") 
write.csv(tnfaSignificant, "resFiles/tnfa-Significant_eGenes_responseEQTL.csv", row.names = F)
write.table(tnfaSignificant[, c(1:4)], "resFiles/tnfa_Significant_eGenes_responseEQTL.txt", row.names = F, col.names = F, quote = F, sep = "\t")


###############################################################################

plot(tnfaUntr$tnfa_slope, tnfaUntr$untr_slope)
plot

ggplot(data = tnfaUntr, aes(y = tnfa_slope, x = untr_slope)) +
  geom_point() +
  labs(title = "Naive eQTL Beta vs TNFa eQTL Beta",
       x = "Naive eQTL Beta",
       y = "TNFa eQTL Beta") +
  theme_minimal() -> tnfaUntrBetaBetaPlot
tnfaUntrBetaBetaPlot
ggsave("tnfaUntrBetaBetaPlot.png", plot = tnfaUntrBetaBetaPlot,
       path = "resPlots", height = 4.73, width = 7.45, units = "in", bg = "white")

ggplot(data = untrTnfa, aes(y = untr_slope, x = tnfa_slope)) +
  geom_point() +
  labs(title = "Naive eQTL Beta vs MAF Plot",
       x = "TNFa eQTL Beta",
       y = "Naive eQTL Beta") +
  theme_minimal() -> untrTnfaBetaBetaPlot
untrTnfaBetaBetaPlot
ggsave("untrTnfaBetaBetaPlot.png", plot = untrTnfaBetaBetaPlot,
       path = "resPlots", height = 4.73, width = 7.45, units = "in", bg = "white")



# Search for assocs
tnfaNom[tnfaNom$tnfa_grp_id == "ENSG00000164308" & tnfaNom$tnfa_rsID == "rs27527", ]
untrNom[untrNom$untr_grp_id == "ENSG00000164308" & untrNom$untr_rsID == "rs27527", ]

##################################################################################
##################################################################################
# sQTL




