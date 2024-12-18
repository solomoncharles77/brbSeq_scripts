
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
tnfaNom <- data.frame(fread("leafcutterSplicing_TNFA_cisEQTL/leafcutterSplicing_TNFA_cisEQTL_nominal_ALL_pval0.05_rsID.txt.gz",
                            select=c(1,2,7,14,16,17,19,24)))
colnames(tnfaNom) <- paste0("tnfa_", colnames(tnfaNom))
gc()

# Import UNTR permute 
untrPem <- data.frame(fread("leafcutterSplicing_UNTR_cisEQTL/leafcutterSplicing_UNTR_cisEQTL_permute_ALL.txt.gz"))
colnames(untrPem) <- paste0("untr_", pemIDs)
untrPem <- untrPem[untrPem$untr_adj_emp_pval < 0.05, ]

tnfaUntr <- merge(tnfaNom, untrPem, by.x = c("tnfa_phe_id", "tnfa_var_id"), by.y = c("untr_phe_id", "untr_var_id"))

tnfaUntr$zScore <- (tnfaUntr$untr_slope - tnfaUntr$tnfa_slope) / sqrt(tnfaUntr$untr_slope_se^2 + tnfaUntr$tnfa_slope_se^2)
#tnfaUntr$zScorePval <- pnorm(tnfaUntr$zScore, mean = 0, sd = 1, lower.tail = TRUE)
tnfaUntr$zScorePval <- 2 * pnorm(-abs(tnfaUntr$zScore))
tnfaUntr2 <- tnfaUntr[tnfaUntr$zScorePval < 0.05, ]
write.csv(tnfaUntr, "resFiles/response_TNFA_vs_UNTR_sQTL.csv", row.names = F)


###########################################################################################

untrNom <- data.frame(fread("geneExprUNTR_brbSeq_cisEQTL/geneExprUNTR_brbSeq_cisEQTL_nominal_ALL_pval0.05_rsID.txt.gz",
                            select=c(1,2,7,14,16,17,19,24)))
colnames(untrNom) <- paste0("untr_", colnames(untrNom))
gc()

tnfaPem <- data.frame(fread("geneExprTNFA_brbSeq_cisEQTL/geneExprTNFA_brbSeq_cisEQTL_permute_ALL.txt.gz"))
colnames(tnfaPem) <- paste0("tnfa_", pemIDs)
tnfaPem <- tnfaPem[tnfaPem$tnfa_adj_emp_pval < 0.05, ]

untrTnfa <- merge(untrNom, tnfaPem, by.x = c("untr_grp_id", "untr_var_id"), by.y = c("tnfa_phe_id", "tnfa_var_id"))

untrTnfa$zScore <- (untrTnfa$untr_slope - untrTnfa$tnfa_slope) / sqrt(untrTnfa$untr_slope_se^2 + untrTnfa$tnfa_slope_se^2)
#untrTnfa$zScorePval <- pnorm(untrTnfa$zScore, mean = 0, sd = 1, lower.tail = TRUE)
untrTnfa$zScorePval <- 2 * pnorm(-abs(untrTnfa$zScore))
untrTnfa2 <- untrTnfa[untrTnfa$zScorePval < 0.05, ]
write.csv(untrTnfa, "resFiles/response_UNTR_vs_TNFA_sQTL.csv", row.names = F)


# Detecting reQTLs by eQTL β-comparison --------------------------------------------------------------

untrTnfa <- read.csv("resFiles/response_UNTR_vs_TNFA_sQTL.csv")
tnfaUntr <- read.csv("resFiles/response_TNFA_vs_UNTR_sQTL.csv")

###############################################################################

plot(tnfaUntr$tnfa_slope, tnfaUntr$untr_slope)
plot

ggplot(data = tnfa, aes(y = tnfa_slope, x = maf)) +
  geom_point() +
  labs(title = "TNFa eQTL Beta vs MAF Plot",
       x = "Minor Allele Frequency",
       y = "Effect Size") +
  theme_minimal() -> tnfaEqtlMafBetaPlot
tnfaEqtlMafBetaPlot
ggsave("tnfaEqtlMafBetaPlot.png", plot = tnfaEqtlMafBetaPlot,
       path = "resPlots", height = 4.73, width = 7.45, units = "in", bg = "white")

ggplot(data = untr, aes(y = untr_slope, x = maf)) +
  geom_point() +
  labs(title = "Naive eQTL Beta vs MAF Plot",
       x = "Minor Allele Frequency",
       y = "Effect Size") +
  theme_minimal() -> naiveEqtlMafBetaPlot
naiveEqtlMafBetaPlot
ggsave("naiveEqtlMafBetaPlot.png", plot = naiveEqtlMafBetaPlot,
       path = "resPlots", height = 4.73, width = 7.45, units = "in", bg = "white")



##################################################################################
##################################################################################
# sQTL




