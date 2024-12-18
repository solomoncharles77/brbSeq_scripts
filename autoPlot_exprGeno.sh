#!/bin/bash

#SBATCH --job-name=plotExprGeno
#SBATCH --account=vasccell
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=40G
#SBATCH --time=23:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cs806@leicester.ac.uk
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --export=NONE

cd $SLURM_SUBMIT_DIR


module load R

# Plot Expr vs Geno with rsID
# TNFA
Rscript brbSeq_scripts/autoPlot_exprGeno.R -c ENSG00000204802 -r rs11263339 -p geneExprTNFA -e phenoFiles/normGeneExpr_TNFA_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
Rscript brbSeq_scripts/autoPlot_exprGeno.R -c ENSG00000254701 -r rs572787392 -p geneExprTNFA -e phenoFiles/normGeneExpr_TNFA_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed


# NAIVE
Rscript brbSeq_scripts/autoPlot_exprGeno.R -c ENSG00000204802 -r rs11263339 -p geneExprNAIVE -e phenoFiles/normGeneExpr_UNTR_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
Rscript brbSeq_scripts/autoPlot_exprGeno.R -c ENSG00000254701 -r rs572787392 -p geneExprNAIVE -e phenoFiles/normGeneExpr_UNTR_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed


# With coord ID
Rscript brbSeq_scripts/autoPlot_exprGeno2.R -i ENSG00000003147 -c 7:8396313:G:T -r rs6966589 -p geneExprTNFA -e phenoFiles/normGeneExpr_TNFA_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
Rscript brbSeq_scripts/autoPlot_exprGeno2.R -i ENSG00000003147 -c 7:8396313:G:T -r rs6966589 -p geneExprNAIVE -e phenoFiles/normGeneExpr_UNTR_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed

Rscript brbSeq_scripts/autoPlot_exprGeno2.R -i ENSG00000106367 -c 7:101451705:A:T -r rs6975978 -p geneExprTNFA -e phenoFiles/geneExpr_TMM_IVTnorm_TNFA_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
Rscript brbSeq_scripts/autoPlot_exprGeno2.R -i ENSG00000106367 -c 7:101451705:A:T -r rs6975978 -p geneExprNAIVE -e phenoFiles/geneExpr_TMM_IVTnorm_UNTR_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed


# For Figures
# AP1S1
Rscript brbSeq_scripts/autoPlot_exprGeno2_TNFA.R -i ENSG00000106367 -c 7:101451705:A:T -r rs6975978 -p geneExprTNFA -e phenoFiles/geneExpr_TMM_IVTnorm_TNFA_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
Rscript brbSeq_scripts/autoPlot_exprGeno2_UNTR.R -i ENSG00000106367 -c 7:101451705:A:T -r rs6975978 -p geneExprNAIVE -e phenoFiles/geneExpr_TMM_IVTnorm_UNTR_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed


# EPC2
Rscript brbSeq_scripts/autoPlot_exprGeno2_TNFA.R -i ENSG00000135999 -c 7:101451705:A:T -r rs6975978 -p geneExprTNFA -e phenoFiles/geneExpr_TMM_IVTnorm_TNFA_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
Rscript brbSeq_scripts/autoPlot_exprGeno2_UNTR.R -i ENSG00000135999 -c 7:101451705:A:T -r rs6975978 -p geneExprNAIVE -e phenoFiles/geneExpr_TMM_IVTnorm_UNTR_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed


# CIZ1
Rscript brbSeq_scripts/autoPlot_exprGeno2_TNFA.R -i ENSG00000148337 -c 9:128178543:C:T -r rs7847709 -p geneExprTNFA -e phenoFiles/geneExpr_TMM_IVTnorm_TNFA_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
Rscript brbSeq_scripts/autoPlot_exprGeno2_UNTR.R -i ENSG00000148337 -c 9:128178543:C:T -r rs7847709 -p geneExprNAIVE -e phenoFiles/geneExpr_TMM_IVTnorm_UNTR_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed

# URI1 untreated
Rscript brbSeq_scripts/autoPlot_exprGeno2_TNFA.R -i ENSG00000105176 -c 19:29875895:C:A -r rs11670202 -p geneExprTNFA -e phenoFiles/geneExpr_TMM_IVTnorm_TNFA_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
Rscript brbSeq_scripts/autoPlot_exprGeno2_UNTR.R -i ENSG00000105176 -c 19:29875895:C:A -r rs11670202 -p geneExprNAIVE -e phenoFiles/geneExpr_TMM_IVTnorm_UNTR_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed

# UR1 tnfa
Rscript brbSeq_scripts/autoPlot_exprGeno2_TNFA.R -i ENSG00000105176 -c 19:29859232:G:A -r rs957165 -p geneExprTNFA -e phenoFiles/geneExpr_TMM_IVTnorm_TNFA_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
Rscript brbSeq_scripts/autoPlot_exprGeno2_UNTR.R -i ENSG00000105176 -c 19:29859232:G:A -r rs957165 -p geneExprNAIVE -e phenoFiles/geneExpr_TMM_IVTnorm_UNTR_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed










# Plot Expr vs Geno

while IFS=$'\t' read -r col1 col2 col3 col4
do

  #Rscript brbSeq_scripts/autoPlot_exprGeno2.R -i ${col1} -c ${col3} -r ${col4} -p geneExprTNFA -e phenoFiles/normGeneExpr_TNFA_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
  #Rscript brbSeq_scripts/autoPlot_exprGeno2.R -i ${col1} -c ${col3} -r ${col4} -p geneExprNAIVE -e phenoFiles/normGeneExpr_UNTR_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
  
  Rscript brbSeq_scripts/autoPlot_exprGeno2.R -i ${col1} -c ${col3} -r ${col4} -p geneExprTNFA -e phenoFiles/geneExpr_TMM_IVTnorm_TNFA_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
  Rscript brbSeq_scripts/autoPlot_exprGeno2.R -i ${col1} -c ${col3} -r ${col4} -p geneExprNAIVE -e phenoFiles/geneExpr_TMM_IVTnorm_UNTR_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed

done < resFiles/sharedSimilar_assocs2.txt	


while IFS=$'\t' read -r col1 col2 col3 col4
do

  # Rscript brbSeq_scripts/autoPlot_exprGeno2.R -i ${col1} -c ${col3} -r ${col4} -p geneExprTNFA -e phenoFiles/normGeneExpr_TNFA_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
  # Rscript brbSeq_scripts/autoPlot_exprGeno2.R -i ${col1} -c ${col3} -r ${col4} -p geneExprNAIVE -e phenoFiles/normGeneExpr_UNTR_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
  
  Rscript brbSeq_scripts/autoPlot_exprGeno2.R -i ${col1} -c ${col3} -r ${col4} -p geneExprTNFA -e phenoFiles/geneExpr_TMM_IVTnorm_TNFA_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
  Rscript brbSeq_scripts/autoPlot_exprGeno2.R -i ${col1} -c ${col3} -r ${col4} -p geneExprNAIVE -e phenoFiles/geneExpr_TMM_IVTnorm_UNTR_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed

done < resFiles/sharedContra_assocs.txt	



while IFS=$'\t' read -r col1 col2 col3 col4
do

  Rscript brbSeq_scripts/autoPlot_exprGeno2.R -i ${col1} -c ${col3} -r ${col4} -p geneExprTNFA -e phenoFiles/normGeneExpr_TNFA_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
  Rscript brbSeq_scripts/autoPlot_exprGeno2.R -i ${col1} -c ${col3} -r ${col4} -p geneExprNAIVE -e phenoFiles/normGeneExpr_UNTR_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed

done < resFiles/tnfa_specific_eGenes_responseEQTL.txt	


while IFS=$'\t' read -r col1 col2 col3 col4
do

  Rscript brbSeq_scripts/autoPlot_exprGeno2.R -i ${col1} -c ${col3} -r ${col4} -p geneExprTNFA -e phenoFiles/normGeneExpr_TNFA_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
  Rscript brbSeq_scripts/autoPlot_exprGeno2.R -i ${col1} -c ${col3} -r ${col4} -p geneExprNAIVE -e phenoFiles/normGeneExpr_UNTR_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed

done < resFiles/naive-specific_eGenes.txt	



# Plot Splicing vs Geno
# # TNFA
# Rscript brbSeq_scripts/autoPlot_exprGeno.R -c ENSG00000050820 -r rs11648176 -p spliceTNFA -e phenoFiles/LeafcutterSplicing_TNFA_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed
# 
# 
# # NAIVE
# Rscript brbSeq_scripts/autoPlot_exprGeno.R -c ENSG00000050820 -r rs11648176 -p spliceNAIVE -e phenoFiles/LeafcutterSplicing_UNTR_qtlTools_Ready.bed.gz -g genoFiles/brbSeqGeno_QCed


# Check rsID
#grep -wF "9:61667843:G:C" /scratch/vasccell/cs806/colocalization/dbSNP/hg38.snp151_All.bed
#grep -wF "5:69560101:G:C" /scratch/vasccell/cs806/colocalization/dbSNP/hg38.snp151_All.bed
