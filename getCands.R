library(data.table)

# Import data -------------------------------------------------------------
naiveCondSpec <- read.csv("resFiles/naive_condition_eQTL.csv")

# Compare TNFA nominal with UNTR perm -------------------------------------
tnfaNom <- data.frame(fread("geneExprTNFA_brbSeq_cisEQTL/geneExprTNFA_brbSeq_cisEQTL_nominal_ALL_pval0.05_rsID.txt.gz",
                            select=c(1,2,14,16,17,19,24,25)))
colnames(tnfaNom) <- paste0("tnfa_", colnames(tnfaNom))
gc()

# tnfaNom2 <- tnfaNom[tnfaNom$tnfa_grp_id %in% naiveCondSpec$geneID & tnfaNom$tnfa_var_id %in% naiveCondSpec$var_id, ]
# tnfaNom2 <- tnfaNom2[tnfaNom2$tnfa_nom_pval > 0.05, ]

tnfaNom2 <- merge(tnfaNom, naiveCondSpec, by.x = c("tnfa_grp_id", "tnfa_var_id"), by.y = c("geneID", "var_id"))

ncsCands <- tnfaNom2[!duplicated(tnfaNom2$tnfa_grp_id), c(1,8,2,6)]

write.table(ncsCands, "resFiles/ncsCandsPlot.txt", sep = "\t", quote = F, row.names = F, col.names = F)


#########################################################################################################
tnfaRespSpec <- read.csv("resFiles/tnfa_response_eQTL.csv")

untrNom <- data.frame(fread("geneExprUNTR_brbSeq_cisEQTL/geneExprUNTR_brbSeq_cisEQTL_nominal_ALL_pval0.05_rsID.txt.gz",
                            select=c(1,2,14,16,17,19,24,25)))
colnames(untrNom) <- paste0("untr_", colnames(untrNom))
gc()

untrNom2 <- untrNom[untrNom$untr_grp_id %in% naiveCondSpec$geneID & untrNom$untr_var_id %in% tnfaRespSpec$var_id, ]
untrNom2 <- untrNom2[untrNom2$untr_nom_pval > 0.05, ]
trsCands <- untrNom2[!duplicated(untrNom2$untr_grp_id), c(1,8)]

write.table(trsCands, "resFiles/trsCands.txt", sep = "\t", quote = F, row.names = F, col.names = F)
