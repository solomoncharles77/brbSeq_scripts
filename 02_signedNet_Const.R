
# Load libraries and make settings-----------------------------------------
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Import gene expression and phenotype data -------------------------------
load("wgcnaFiles/huvecTNFA_WGCNA_dataInput.RData")

# Network construction ----------------------------------------------------
signed_allExpr_Net = blockwiseModules(datExpr, power = 3, networkType = "signed", maxBlockSize = 20000,
                                      TOMType = "unsigned", minModuleSize = 30, corType = "bicor",
                                      mergeCutHeight = 0.25,
                                      numericLabels = TRUE, pamRespectsDendro = FALSE,
                                      saveTOMs = TRUE,
                                      saveTOMFileBase = "huvecTNFA_WGCNA_signedNet_TOM",
                                      verbose = 5)

save(signed_allExpr_Net, file = "wgcnaFiles/huvecTNFA_WGCNA_signedNet.RData")

rm(list = ls())
gc()
##############################################################

# Import gene expression and phenotype data -------------------------------
load("wgcnaFiles/huvecUNTR_WGCNA_dataInput.RData")


# Network construction ----------------------------------------------------
signed_allExpr_Net = blockwiseModules(datExpr, power = 3, networkType = "signed", maxBlockSize = 30000,
                                      TOMType = "unsigned", minModuleSize = 30, corType = "bicor",
                                      mergeCutHeight = 0.25,
                                      numericLabels = TRUE, pamRespectsDendro = FALSE,
                                      saveTOMs = TRUE,
                                      saveTOMFileBase = "VSMC_AllFemale_Net_TOM",
                                      verbose = 5)

save(signed_allExpr_Net, file = "wgcnaFiles/huvecUNTR_WGCNA_signedNet.RData")
rm(list = ls())
gc()

# #####################################################################################
# # Import gene expression and phenotype data -------------------------------
# load("VSMC_Consensus-dataInput.RData")
# 
# 
# net = blockwiseConsensusModules(
#   multiExpr, power = 9, minModuleSize = 30, deepSplit = 2,
#   pamRespectsDendro = FALSE, maxBlockSize = 25000, corType = "bicor",
#   mergeCutHeight = 0.25, numericLabels = TRUE,
#   minKMEtoStay = 0,
#   saveTOMs = TRUE, verbose = 5)
# 
# consMEs = net$multiMEs;
# moduleLabels = net$colors;
# # Convert the numeric labels to color labels
# moduleColors = labels2colors(moduleLabels)
# consTree = net$dendrograms[[1]]; 
# 
# 
# sizeGrWindow(8,6);
# pdf(file = "Plots/VSMC_ConsensusDendrogram-auto.pdf", wi = 8, he = 6)
# plotDendroAndColors(consTree, moduleColors,
#                     "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05,
#                     main = "Consensus gene dendrogram and module colors")
# 
# dev.off()
# 
# 
# save(net, consMEs, moduleLabels, moduleColors, consTree, file = "VSMC_Consensus-NetworkConstruction-auto.RData")
