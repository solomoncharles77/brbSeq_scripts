# Count the number of significant eQTLs from fastQTL 
# permutation pass. 
# Boxiang Liu
# 2018-02-14

library(data.table)
library(foreach)
in_fn = "geneExprTNFA_brbSeq_cisEQTL/geneExprTNFA_brbSeq_cisEQTL_permute_ALL.txt.gz"
out_dir = 'resFiles'

#if (!dir.exists(out_dir)) {dir.create(out_dir)}

read_eqtl = function(in_fn){
  eqtl = fread(in_fn,select=c(1,6,9,22),col.names=c('grp_id','phe_id','dist_phe_var','pval'))
  return(eqtl)
}

adjust_pvalue = function(pval){
  padj = p.adjust(pval,method='fdr')
  return(padj)
}


count_sig_eqtl = function(padj,levels = c(0.05,0.01,0.001)){
  sig_eqtl = foreach (level = levels,.combine='rbind')%do%{
    n_sig = sum(padj < level, na.rm = TRUE)
    data.table(
      level = level,
      sig = n_sig
    )
  }
  return(sig_eqtl)
}

save_result = function(output, out_fn){
  fwrite(
    output,
    out_fn,
    sep = '\t'
  )
}
eqtl = read_eqtl(in_fn)
eqtl$padj = adjust_pvalue(eqtl$pval)
sig_eqtl = count_sig_eqtl(eqtl$padj)
sig_eqtl
save_result(sig_eqtl,sprintf('%s/sig_eqtl.tsv',out_dir))