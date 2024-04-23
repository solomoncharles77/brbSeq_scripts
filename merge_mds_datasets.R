# Merge read count data for MDS plot
# Charles Solomon
# 05/05/2020


# Module load libraries ---------------------------------------------------
library(data.table)
library(biomaRt)

# Import data sets --------------------------------------------------------
# ourHUVEC
ourHUVEC <- data.frame(fread("/scratch/cellfunc/cs806/huvecMDS/dataFiles/ourHUVEC_Gene_Expr_Counts.txt"))

# Get gene information with Biomart
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
geneInfo <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "start_position","end_position"),
                  filters="ensembl_gene_id", values=ourHUVEC$GeneID, mart=human)
geneInfo$size <- geneInfo$end_position - geneInfo$start_position

ourHUVEC <- merge(geneInfo[, c(1,2)], ourHUVEC, by.x = "ensembl_gene_id", by.y = "GeneID")
ourHUVEC <- ourHUVEC[, -c(1,3)]
colnames(ourHUVEC)[1] <- "GeneID"
ourHUVEC <- ourHUVEC[!(duplicated(ourHUVEC$GeneID) | duplicated(ourHUVEC$GeneID, fromLast = TRUE)), ] # get rid of duplicated gene names
ourHUVEC[1:10, 1:10]

# BRBseq
untrBRBseq <- data.frame(fread("exprFiles/untreated_pilot_BRBseq.txt"))
untrBRBseq <- merge(geneInfo[, c(1,2)], untrBRBseq, by.x = "ensembl_gene_id", by.y = "GeneID")
untrBRBseq <- untrBRBseq[, -c(1)]
colnames(untrBRBseq)[1] <- "GeneID"
colnames(untrBRBseq)[-1] <- paste0(colnames(untrBRBseq)[-1], "_untr")
untrBRBseq <- untrBRBseq[!(duplicated(untrBRBseq$GeneID) | duplicated(untrBRBseq$GeneID, fromLast = TRUE)), ] # get rid of duplicated gene names
untrBRBseq[1:10, 1:10]


# tnfaBRBseq
tnfaBRBseq <- data.frame(fread("exprFiles/tnfa_pilot_BRBseq.txt"))
tnfaBRBseq <- merge(geneInfo[, c(1,2)], tnfaBRBseq, by.x = "ensembl_gene_id", by.y = "GeneID")
tnfaBRBseq <- tnfaBRBseq[, -c(1)]
colnames(tnfaBRBseq)[1] <- "GeneID"
colnames(tnfaBRBseq)[-1] <- paste0(colnames(tnfaBRBseq)[-1], "_tnfa")
tnfaBRBseq <- tnfaBRBseq[!(duplicated(tnfaBRBseq$GeneID) | duplicated(tnfaBRBseq$GeneID, fromLast = TRUE)), ] # get rid of duplicated gene names
tnfaBRBseq[1:10, 1:10]


# hcaec
hcaec <- data.frame(fread("/scratch/cellfunc/cs806/huvecMDS/dataFiles/HCAEC_expression.txt"))
colnames(hcaec)[1] <- "GeneID"
hcaec <- hcaec[!(duplicated(hcaec$GeneID) | duplicated(hcaec$GeneID, fromLast = TRUE)), ] # get rid of duplicated gene names
hcaec[1:10, 1:10]

# haec
haec <- data.frame(fread("/scratch/cellfunc/cs806/huvecMDS/dataFiles/HAEC_expression.txt"))
colnames(haec)[1] <- "GeneID"
haec <- haec[!(duplicated(haec$GeneID) | duplicated(haec$GeneID, fromLast = TRUE)), ] # get rid of duplicated gene names
haec[1:10, 1:10]

# hcec
hcec <- data.frame(fread("/scratch/cellfunc/cs806/huvecMDS/dataFiles/HCEC_expression.txt"))
colnames(hcec)[1] <- "GeneID"
hcec <- hcec[!(duplicated(hcec$GeneID) | duplicated(hcec$GeneID, fromLast = TRUE)), ] # get rid of duplicated gene names
hcec[1:10, 1:10]

# beta
beta <- data.frame(fread("/scratch/cellfunc/cs806/huvecMDS/dataFiles/BETA_expression.txt"))
colnames(beta)[1] <- "GeneID"
beta <- beta[!(duplicated(beta$GeneID) | duplicated(beta$GeneID, fromLast = TRUE)), ] # get rid of duplicated gene names
beta[1:10, 1:10]

# fibroblast
fibroblast <- data.frame(fread("/scratch/cellfunc/cs806/huvecMDS/dataFiles/FIBROBLASTS_expression.txt"))
colnames(fibroblast)[1] <- "GeneID"
fibroblast <- fibroblast[!(duplicated(fibroblast$GeneID) | duplicated(fibroblast$GeneID, fromLast = TRUE)), ] # get rid of duplicated gene names
fibroblast[1:10, 1:10]

# neutrophil
neutrophil <- data.frame(fread("/scratch/cellfunc/cs806/huvecMDS/dataFiles/NEUTROPHIL_expression.txt"))
colnames(neutrophil)[1] <- "GeneID"
neutrophil <- neutrophil[!(duplicated(neutrophil$GeneID) | duplicated(neutrophil$GeneID, fromLast = TRUE)), ] # get rid of duplicated gene names
neutrophil[1:10, 1:10]

# macrophage
macrophage <- data.frame(fread("/scratch/cellfunc/cs806/huvecMDS/dataFiles/MACROPHAGE_expression.txt"))
colnames(macrophage)[1] <- "GeneID"
macrophage <- macrophage[!(duplicated(macrophage$GeneID) | duplicated(macrophage$GeneID, fromLast = TRUE)), ] # get rid of duplicated gene names
macrophage[1:10, 1:10]

# huvec
huvec <- data.frame(fread("/scratch/cellfunc/cs806/huvecMDS/dataFiles/HUVEC_expression.txt"))
colnames(huvec)[1] <- "GeneID"
huvec <- huvec[!(duplicated(huvec$GeneID) | duplicated(huvec$GeneID, fromLast = TRUE)), ] # get rid of duplicated gene names
huvec[1:10, 1:10]


# merge  and export all data sets
# df_list <- list(ourHUVEC, hcaec, haec, hcec, hcasmc)
# df_list <- list(ourHUVEC, hcaec, haec, hcec, hcasmc, beta, fibroblast, macrophage, neutrophil)
# df_list <- list(ourHUVEC, hcaec, haec, hcec, hcasmc, beta, fibroblast, macrophage, neutrophil, huvec)
df_list <- list(ourHUVEC, hcaec, haec, hcec, fibroblast, macrophage, neutrophil, huvec, untrBRBseq, tnfaBRBseq)
ecDF <- Reduce(function(x, y) merge(x, y, all=FALSE), df_list)
ecDF[1:10, 1:10]
fwrite(ecDF, "mdsFiles/merged_expression_count_4_MDS_plot.txt")


# Import, create, merge and export attribute files ----------------------------------------
hcaeAtt <- read.delim("/scratch/cellfunc/cs806/huvecMDS/dataFiles/HCAEC_attributes.txt")
haecAtt <- read.delim("/scratch/cellfunc/cs806/huvecMDS/dataFiles/HAEC_attributes.txt")
hcecAtt <- read.delim("/scratch/cellfunc/cs806/huvecMDS/dataFiles/HCEC_attributes.txt")

ourHUVEC_Att <- data.frame(sampID = colnames(ourHUVEC)[-c(1)], sampName = "ourHUVEC",
                           sampTitle = "ourHUVEC", sampChar = "ourHUVEC", alias = "ourHUVEC")
untrBRBseq_Att <- data.frame(sampID = colnames(untrBRBseq)[-c(1)], sampName = "untrBRBseq",
                             sampTitle = "untrBRBseq", sampChar = "untrBRBseq", alias = "untrBRBseq")
tnfaBRBseq_Att <- data.frame(sampID = colnames(tnfaBRBseq)[-c(1)], sampName = "tnfaBRBseq",
                           sampTitle = "tnfaBRBseq", sampChar = "tnfaBRBseq", alias = "tnfaBRBseq")

# hcasmcAtt <- data.frame(sampID = colnames(hcasmc)[-c(1)], sampName = "HCASMC",
#                         sampTitle = "HCASMC", sampChar = "HCASMC", alias = "HCASMC")

betaAtt <- read.delim("/scratch/cellfunc/cs806/huvecMDS/dataFiles/BETA_attributes.txt")
betaAtt$sampContact <- NULL
fibroblastAtt <- read.delim("/scratch/cellfunc/cs806/huvecMDS/dataFiles/FIBROBLASTS_attributes.txt")
fibroblastAtt$sampContact <- NULL
neutrophilAtt <- read.delim("/scratch/cellfunc/cs806/huvecMDS/dataFiles/NEUTROPHIL_attributes.txt")
neutrophilAtt$sampContact <- NULL
macrophageAtt <- read.delim("/scratch/cellfunc/cs806/huvecMDS/dataFiles/MACROPHAGE_attributes.txt")
macrophageAtt$sampContact <- NULL
huvecAtt <- read.delim("/scratch/cellfunc/cs806/huvecMDS/dataFiles/HUVEC_attributes.txt")
huvecAtt$sampContact <- NULL

# df_list2 <- list(ourHUVEC_Att, hcaeAtt, haecAtt, hcecAtt, hcasmcAtt)

# df_list2 <- list(ourHUVEC_Att, hcaeAtt, haecAtt, hcecAtt, hcasmcAtt,
#                  betaAtt, fibroblastAtt,  macrophageAtt, neutrophilAtt)

# df_list2 <- list(ourHUVEC_Att, hcaeAtt, haecAtt, hcecAtt, hcasmcAtt,
#                  betaAtt, fibroblastAtt,  macrophageAtt, neutrophilAtt, huvecAtt)

df_list2 <- list(ourHUVEC_Att,  hcaeAtt, haecAtt, hcecAtt,
                 fibroblastAtt,  macrophageAtt, neutrophilAtt, huvecAtt, untrBRBseq_Att, tnfaBRBseq_Att)

ecAtt <- rbindlist(df_list2)
ecAtt[1:10, 1:5]
fwrite(ecAtt, "mdsFiles/merged_expression_count_sample_attributes_4_MDS_plot.txt")
