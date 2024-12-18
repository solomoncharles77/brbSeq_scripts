library(data.table)
library(tidyverse)

# Write function
searchAssoc <- function(qtl, geneName, rsID){
  df <- qtl[qtl$external_gene_name == geneName, ]
  df <- df[df$rsID == rsID, ]
  return(df)
  
}


# Compare TNFA nominal with UNTR perm -------------------------------------
tnfaNom <- data.frame(fread("geneExprTNFA_brbSeq_cisEQTL/geneExprTNFA_brbSeq_cisEQTL_nominal_ALL_pval0.05_rsID.txt.gz"))
gc()

untrNom <- data.frame(fread("geneExprUNTR_brbSeq_cisEQTL/geneExprUNTR_brbSeq_cisEQTL_nominal_ALL_pval0.05_rsID.txt.gz"))
gc()

gn <- "EPC2"
rs <- "rs10754978"

adnp2_tnfa <- searchAssoc(tnfaNom, gn, rs)
adnp2_untr <- searchAssoc(untrNom, gn, rs)
adnp2_tnfa
adnp2_untr


tnfaExpr <- data.frame(fread("phenoFiles/geneExpr_TMM_IVTnorm_TNFA_qtlTools_Ready.bed.gz"))
untrExpr <- data.frame(fread("phenoFiles/geneExpr_TMM_IVTnorm_UNTR_qtlTools_Ready.bed.gz"))

tnfaExpr[1:10, 1:10]
tnfaExpr[tnfaExpr$gene == "ENSG00000101544", ]

untrExpr[untrExpr$gene == "ENSG00000101544", ]

# EPC2
Rscript brbSeq_scripts/autoPlot_exprGeno2_TNFA.R -i ENSG00000135999 -c 2:149158515:A:G -r rs10754978 -p geneExprTNFA -e phenoFiles/geneExpr_TMM_IVTnorm_TNFA_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
Rscript brbSeq_scripts/autoPlot_exprGeno2_UNTR.R -i ENSG00000135999 -c 2:149158515:A:G -r rs10754978 -p geneExprNAIVE -e phenoFiles/geneExpr_TMM_IVTnorm_UNTR_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed

