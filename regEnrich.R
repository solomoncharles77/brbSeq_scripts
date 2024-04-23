
library(RegEnrich)
library(data.table)

expr <- data.frame(fread("HUVEC_BRBSeq_Pilot_Study/AMP0157_results/count_matrix/L101123LU01_01.read.counts.sampleIDs.txt"))
row.names(expr) <- expr$Gene_id
expr$Gene_id <- NULL
expr <- expr[ rowSums(expr) > ncol(expr), ]

sampInfo <- data.frame(sampID = colnames(expr),
                       treatment = sub(".*_", "", colnames(expr)),
                       source_name = "HUVEC",
                       organism = "Homo sapiens",
                       molecule = "total RNA")

data(TFs)
designMatrix = ~ treatment 
contrast = c(0, 1)

object = RegenrichSet(expr = expr, # expression data (matrix)
                      colData = sampInfo, # sample information (data frame)
                      reg = TFs$TF, # regulators
                      method = "Wald_DESeq2", # differentila expression analysis method
                      minMeanExpr = ncol(expr),
                      design = designMatrix, # desing model matrix
                      contrast = contrast, # contrast
                      networkConstruction = "COEN", # network inference method
                      enrichTest = "FET") # enrichment analysis method  


# Differential expression analysis
object = regenrich_diffExpr(object)
# Regulator-target network inference
object = regenrich_network(object)
# Enrichment analysis
object = regenrich_enrich(object)
# Regulator scoring and ranking
object = regenrich_rankScore(object)

# Obtaining results
res = results_score(object)





# Perform 4 steps in RegEnrich analysis
object = regenrich_diffExpr(object) %>% 
  regenrich_network() %>% 
  regenrich_enrich() %>% 
  regenrich_rankScore()

res = results_score(object)