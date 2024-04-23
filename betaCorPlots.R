
library(data.table)
library(tidyverse)
library(openxlsx)
library(gridExtra)


# Import datasets --------------------------------------------------------
huvec50 <- data.frame(fread("huvec50_cisEQTL_2pc_3pf_Results/huvec50_cisEQTL_2pc_3pf.csv"))
huvec50 <- huvec50[huvec50$pvalue < 0.05, ]
tnfa <- data.frame(fread("tnfa_pilot_BRBseq_cisEQTL_2pc_3pf_Results/tnfa_pilot_BRBseq_cisEQTL_2pc_3pf.csv"))
tnfa <- tnfa[tnfa$pvalue < 0.05, ]
untr <- data.frame(fread("untreated_pilot_BRBseq_cisEQTL_2pc_3pf_Results/untreated_pilot_BRBseq_cisEQTL_2pc_3pf.csv"))
untr <- untr[untr$pvalue < 0.05, ]


# Compare bulk HUVEC and Control BRBSeq HUVEC -----------------------------
huvec50_untr <- merge(huvec50, untr, by = c("gene", "snps"), suffixes = c("_huvec50", "_brbSeqUntr"))
huvec50_untr$chr <- sapply(huvec50_untr$snps, function(x){unlist(strsplit(x, split = ":|_"))[1]})
huvec50_untr <- split(huvec50_untr, huvec50_untr$chr)
write.xlsx(huvec50_untr, file = "corPlots/huvec50_untreated_BRBseq_pilot_compare_pvalue_0.05.xlsx")

huvec50_untrPlots <- lapply(names(huvec50_untr),function(chr) {
  df <- huvec50_untr[[chr]]
  df <- df[, c(6,10)]
  corr <- cor(df, use="complete.obs", method="spearman")
  corr <- format(corr, digits=3, nsmall=3)
  
  p<-ggplot(data=df, aes(x=beta_huvec50, y=beta_brbSeqUntr)) +
    geom_point()
  p <- p + stat_smooth(method="lm") + 
    labs(title = chr) +
    theme_linedraw(base_size=10) +
    annotate("text", x=-1.3, y=1.3, label=paste0("r = ", corr[1,2]))
})

huvec50_untr_betaCorPlot <- marrangeGrob(huvec50_untrPlots, nrow=3, ncol=3, top = "")
ggsave(paste0("huvec50BETA_Vs_untreatedBRBseqBETA_Corplots_pvalue_0.05.pdf"), huvec50_untr_betaCorPlot,
       path = "corPlots",
       width = 20, height = 20, units = "cm",
       limitsize = TRUE)

# Compare bulk HUVEC and TNFA BRBSeq HUVEC -----------------------------
huvec50_tnfa <- merge(huvec50, tnfa, by = c("gene", "snps"), suffixes = c("_huvec50", "_brbSeqTNFA"))
huvec50_tnfa$chr <- sapply(huvec50_tnfa$snps, function(x){unlist(strsplit(x, split = ":|_"))[1]})
huvec50_tnfa <- split(huvec50_tnfa, huvec50_tnfa$chr)
write.xlsx(huvec50_tnfa, file = "corPlots/huvec50_TNFA_BRBseq_pilot_compare_pvalue_0.05.xlsx")

huvec50_tnfaPlots <- lapply(names(huvec50_tnfa),function(chr) {
  df <- huvec50_tnfa[[chr]]
  df <- df[, c(6,10)]
  corr <- cor(df, use="complete.obs", method="spearman")
  corr <- format(corr, digits=3, nsmall=3)
  
  p<-ggplot(data=df, aes(x=beta_huvec50, y=beta_brbSeqTNFA)) +
    geom_point()
  p <- p + stat_smooth(method="lm") + 
    labs(title = chr) +
    theme_linedraw(base_size=10) +
    annotate("text", x=-1.3, y=1.3, label=paste0("r = ", corr[1,2]))
})

huvec50_tnfa_betaCorPlot <- marrangeGrob(huvec50_tnfaPlots, nrow=3, ncol=3, top = "")
ggsave(paste0("huvec50BETA_Vs_TNFA_BETA_Corplots_pvalue_0.05.pdf"), huvec50_tnfa_betaCorPlot,
       path = "corPlots",
       width = 20, height = 20, units = "cm",
       limitsize = TRUE)


# Compare Control and TNFA treated BRBSeq HUVEC  -----------------------------
untr_tnfa <- merge(untr, tnfa, by = c("gene", "snps"), suffixes = c("_brbSeqUntr", "_brbSeqTNFA"))
untr_tnfa$chr <- sapply(untr_tnfa$snps, function(x){unlist(strsplit(x, split = ":|_"))[1]})
untr_tnfa <- split(untr_tnfa, untr_tnfa$chr)
write.xlsx(untr_tnfa, file = "corPlots/Control_TNFA_BRBseq_pilot_compare_pvalue_0.05.xlsx")


untr_tnfaPlots <- lapply(names(untr_tnfa),function(chr) {
  df <- untr_tnfa[[chr]]
  df <- df[, c(6,10)]
  corr <- cor(df, use="complete.obs", method="spearman")
  corr <- format(corr, digits=3, nsmall=3)
  
  p<-ggplot(data=df, aes(x=beta_brbSeqUntr, y=beta_brbSeqTNFA)) +
    geom_point()
  p <- p + stat_smooth(method="lm") + 
    labs(title = chr) +
    theme_linedraw(base_size=10) +
    annotate("text", x=-1.3, y=1.3, label=paste0("r = ", corr[1,2]))
})

untr_tnfa_betaCorPlot <- marrangeGrob(untr_tnfaPlots, nrow=3, ncol=3, top = "")
ggsave(paste0("Control_BRBseq_Vs_TNFA_BETA_Corplots_pvalue_0.05.pdf"), untr_tnfa_betaCorPlot,
       path = "corPlots",
       width = 20, height = 20, units = "cm",
       limitsize = TRUE)
