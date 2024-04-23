library(data.table)

expr <- data.frame(fread("HUVEC_BRBSeq_Pilot_Study/AMP0157_results/count_matrix/L101123LU01_01.read.counts.sampleIDs.txt"))
fam <- read.table("../genotypeData/rawGeno/Imputed_Data_TOPMed_hg38.fam")

# Get and export genotype data
ss <- sub("_.*", "", colnames(expr)[-1])
ss <- sub("X", "S", ss)
ssFam <- fam[fam$V2 %in% ss, ]
write.table(ssFam, "../genotypeData/rawGeno/brbSeq_Samplist.txt", row.names = F, col.names = F, quote = F)

# Split expr into treated and untreated
row.names(expr) <- expr$Gene_id
expr$Gene_id <- NULL
trt <- c("untr", "TNFa")
exprList <- lapply(setNames(trt, trt), function(x) expr[, grep(x, colnames(expr))])

# Organize and export each subset 
untr <- exprList$untr
colnames(untr) <- sub("_.*", "", colnames(untr))
colnames(untr) <- sub("X", "S", colnames(untr))
untr <- cbind(GeneID = rownames(untr), untr)
fwrite(untr, "exprFiles/untreated_pilot_BRBseq.txt", sep = "\t")

tnfa <- exprList$TNFa
colnames(tnfa) <- sub("_.*", "", colnames(tnfa))
colnames(tnfa) <- sub("X", "S", colnames(tnfa))
tnfa <- cbind(GeneID = rownames(tnfa), tnfa)
fwrite(tnfa, "exprFiles/tnfa_pilot_BRBseq.txt", sep = "\t")

rm(list = ls())
#########################################################################################################################

# # Import and process bulk expr data --------------------------------------------
# expr <- read.table("/scratch/cellfunc/shared/HUVEC_BackUp_Oct22/HUVEC_RNAseq/Counts/salmon.merged.gene_counts.tsv", header = TRUE)
# # Harmonize names and samples in geno and expr ----------------------------
# mapFile <- read.csv("/scratch/cellfunc/shared/HUVEC_BackUp_Oct22/HUVEC_RNAseq/Mappings.csv")
# names(expr) <- mapFile$supplier_name[match(names(expr), mapFile$sanger_sample_id)]
# names(expr)[1:2] <- c("ensembl_gene_id", "external_gene_name")
# rownames(expr) <- expr$ensembl_gene_id
# expr <- expr[, -c(1:3)]
# colnames(expr) <- sub("E", "S", colnames(expr))
# 
# fam <- read.table("/scratch/cellfunc/shared/HUVEC_BackUp_Oct22/HUVEC_genotype/Imputed_Genotypes_Data_HUVEC_Samples_with_RNA-seq_Data.fam")
# ssFam <- fam[fam$V2 %in% colnames(expr), ]
# 
# 
# # Sample a random set of males and females 25 individuals each ------------
# set.seed(123)
# 
# # Sample 25 males
# mFam <- ssFam[ssFam$V5 == 1, , drop = FALSE]
# mFam_sampled <- mFam[sample(nrow(mFam), 25), ]
# 
# # Sample 25 females
# fFam <- ssFam[ssFam$V5 == 2, , drop = FALSE]
# fFam_sampled <- fFam[sample(nrow(fFam), 25), ]
# 
# fam50 <- rbind(mFam_sampled, fFam_sampled)
# write.table(fam50, "/scratch/cellfunc/shared/HUVEC_BackUp_Oct22/HUVEC_genotype/huvec50.txt", row.names = F, col.names = F, quote = F)
# 
# expr50 <- expr[, names(expr) %in% fam50$V2]
# expr50 <- cbind(GeneID = rownames(expr50), expr50)
# fwrite(expr50, "exprFiles/huvec50.txt", sep = "\t")

############################################################################
# export result
ff <- c("untreated_pilot_BRBseq_cisEQTL_2pc_3pf_Results/untreated_pilot_BRBseq_cisEQTL_2pc_3pf_QQPlot.png",
        "untreated_pilot_BRBseq_cisEQTL_2pc_3pf_Results/untreated_pilot_BRBseq_cisEQTL_2pc_3pf_4MR.txt.gz",
        "untreated_pilot_BRBseq_cisEQTL_2pc_3pf_Results/untreated_pilot_BRBseq_cisEQTL_2pc_3pf_betaMaf_summary.csv",
        "tnfa_pilot_BRBseq_cisEQTL_2pc_3pf_Results/tnfa_pilot_BRBseq_cisEQTL_2pc_3pf_QQPlot.png",
        "tnfa_pilot_BRBseq_cisEQTL_2pc_3pf_Results/tnfa_pilot_BRBseq_cisEQTL_2pc_3pf_4MR.txt.gz",
        "tnfa_pilot_BRBseq_cisEQTL_2pc_3pf_Results/tnfa_pilot_BRBseq_cisEQTL_2pc_3pf_betaMaf_summary.csv",
        "huvec50_cisEQTL_2pc_3pf_Results/huvec50_cisEQTL_2pc_3pf_QQPlot.png",
        "huvec50_cisEQTL_2pc_3pf_Results/huvec50_cisEQTL_2pc_3pf_4MR.txt.gz",
        "huvec50_cisEQTL_2pc_3pf_Results/huvec50_cisEQTL_2pc_3pf_betaMaf_summary.csv",
        "mafBetaPlots",
        "corPlots")
zip::zip("pilot_BRBSeq_eQTL_data.zip", files = ff)


cc <- c("corPlots/Control_BRBseq_Vs_TNFA_BETA_Corplots_FDR_0.05.pdf",
        "corPlots/Control_BRBseq_Vs_TNFA_BETA_Corplots_pvalue_0.05.pdf",
        "corPlots/huvec50BETA_Vs_TNFA_BETA_Corplots_FDR_0.05.pdf",
        "corPlots/huvec50BETA_Vs_TNFA_BETA_Corplots_pvalue_0.05.pdf",
        "corPlots/huvec50BETA_Vs_untreatedBRBseqBETA_Corplots_FDR_0.05.pdf",
        "corPlots/huvec50BETA_Vs_untreatedBRBseqBETA_Corplots_pvalue_0.05.pdf")
zip::zip("FDR_pvalue_filtered_Corplots.zip", files = cc)

