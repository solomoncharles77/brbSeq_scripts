# This script is run with R v3.6.1
# Code source https://www.biostars.org/p/162566/
# load libraries

# This script takes one argument -  the prefix name of the plink processed output
aha <- commandArgs(trailingOnly = TRUE)

aha2 <- paste0("exprFiles/", aha, "_geneExpr_qtlReady.txt")


library(peer)
library(data.table)


# Import expression -------------------------------------------------------
expr <- read.table(aha2, row.names = 1,  header = T)

# the PEER tutorial recommends k == 25% of "number of samples"
k = ncol(expr) * .25

# create empty peer model
model = PEER()

# set number of "hidden factors" searched for
PEER_setNk(model,k)

# PEER ask NxG matrix, where N=samples and G=genes
PEER_setPhenoMean(model,as.matrix(t(expr)))

PEER_setAdd_mean(model, TRUE)

# These can be changed to the tolerance you desire
#PEER_setTolerance(model, .001)
#PEER_setVarTolerance(model, 0.0005)

# usually doesn't go over 200 iterations        
#PEER_setNmax_iterations(model, 1000)

# run peer
PEER_update(model)

allPeer = PEER_getX(model)


factors = data.frame(t(PEER_getX(model)))

colnames(factors) <- colnames(expr)
factors <- factors[-1, ]
rownames(factors) <- paste0("peerF", 1:nrow(factors))
factors <- data.frame(Covariates = rownames(factors), factors)

peerFileName <- paste0("covFiles/", aha, "_PeerFactor.txt")
fwrite(factors, file = peerFileName, sep = "\t", quote = F)

saveRDS(allPeer, file = paste0("covFiles/", aha, ".rds"))

