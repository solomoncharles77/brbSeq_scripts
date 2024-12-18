
library(tidyverse)
library(data.table)

qtlHeader <- read.delim("brbSeq_scripts/QTLtools_permute_header.txt", header = F)
colIDs <- qtlHeader$V2
colIDs <- sub(" \\| ve_by_pc1 \\| n_phe_in_grp", "", colIDs)
colIDs <- sub("phe_id \\| ", "", colIDs)

maf <- data.frame(fread("genoFiles/brbSeqGeno_QCed_MAF.txt"))
colnames(maf) <- c("var_id", "chr", "pos", "ref", "alt", "maf")

annot <- data.frame(read.csv("exprFiles/geneAnnot.csv"))

addRSID <- function(file_path, feature_name) {
  sumstat <- data.frame(fread(file_path, header = F, stringsAsFactors = F))
  colnames(sumstat) <- colIDs
  
  coordFileOut <- paste0("../colocalization/summStatsBEDs/", feature_name, "_coordID.txt")
  fwrite(sumstat[, c("var_id", "var_id")], coordFileOut, sep = "\t", quote = F, row.names = F, col.names = F)
  coord <- system(paste0("/home/c/cs806/tsv-utils-v2.1.2_linux-x86_64_ldc2/bin/tsv-join -f ", coordFileOut, " -k 1 -d 2 ../colocalization/dbSNP/dbSNP155.hg38.rsID_CoordID.txt"), intern = T)
  coord <- data.frame(fread(text = coord, header = F))
  colnames(coord) <- c("rsID", "var_id")
  
  sumstat2 <- merge(sumstat, coord, by = "var_id", all.x = TRUE)
  gc()
  sumstat2$variance <- sumstat2$slope_se * 2
  sumstat2 <- merge(sumstat2, maf[, c(1,3:6)], by = "var_id", all.x = TRUE)
  gc()
  sumstat2 <- merge(sumstat2, annot[, -3], by.x = "grp_id", by.y = "ensembl_gene_id", all.x = TRUE)
  gc()
  
  return(sumstat2)
}


tnfaEqtl <- addRSID("geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL/geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL_permute_ALL.txt.gz", "tnfaEqtl")
fwrite(tnfaEqtl, "geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL/geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL_permute_ALL_rsID.txt.gz", sep = "\t")
gc()

untrEqtl <- addRSID("geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL/geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL_permute_ALL.txt.gz", "untrEqtl")
fwrite(untrEqtl, "geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL/geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL_permute_ALL_rsID.txt.gz", sep = "\t")
gc()


######################################################################################################################

