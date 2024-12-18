
library(data.table)
library(locuszoomr)
library(AnnotationHub)

ah <- AnnotationHub()
query(ah, c("EnsDb", "v100", "Homo sapiens"))
#ensDb_v106 <- ah[["AH100643"]]
ensDb_v100 <- ah[["AH79689"]]


qtl1 <- data.frame(fread("geneExprUNTR_brbSeq_cisEQTL/geneExprUNTR_brbSeq_cisEQTL_nominal_ALL_pval0.05_rsID.txt.gz",
                         select = c(7,25,11,12,19,22,23,14,16,17,15),
                         col.names = c("pheID", "GeneName", "chrom", "pos", "rsid",
                                       "other_allele", "effect_allele", "p", "beta", "se", "r2")))

qtl2 <- data.frame(fread("geneExprTNFA_brbSeq_cisEQTL/geneExprTNFA_brbSeq_cisEQTL_nominal_ALL_pval0.05_rsID.txt.gz",
                         select = c(7,25,11,12,19,22,23,14,16,17,15),
                         col.names = c("pheID", "GeneName", "chrom", "pos", "rsid",
                                       "other_allele", "effect_allele", "p", "beta", "se", "r2")))
gc()

cand <- "ENSG00000148337"
candName <- "CIZ1"


candQtlSub2 <- qtl2[qtl2$pheID == cand, ]
candQtlSub2 <- candQtlSub2[order(candQtlSub2$p, candQtlSub2$rsid), ]
topSNP2 <- candQtlSub2[1, 5]

snpListFile2 <- paste0(candName, "_SNPList2.txt")
topSNPFile2 <- paste0(candName, "_topSNP2.txt")
snpLDFile2 <- paste0(candName, "_SNPList_LD2.txt")
write.table(data.frame(candQtlSub2$rsid), file = snpListFile2,  row.names = F, col.names = F, quote = F)
write.table(data.frame(topSNP2), file = topSNPFile2,  row.names = F, col.names = F, quote = F)
system(paste0("plink --bfile /scratch/vasccell/cs806/colocalization/1000Genome/g1000_eur --extract ", snpListFile2, " --ld-window-r2 0.00000001 --r2  --ld-window-kb 2000000 --ld-window 1000000 --ld-snp-list ", topSNPFile2, " --threads 6 --out ", snpLDFile2))
ld2 <- read.table(paste0(snpLDFile2, ".ld"), header = T)

candQtlSub2 <- merge(candQtlSub2, ld2[, c(6,7)], by.x = "rsid", by.y = "SNP_B")
candQtlSub2$logP <- -log10(candQtlSub2$p)
candQtlSub2$ld <- candQtlSub2$R2
loc2 <- locus(data = candQtlSub2, gene = candName, fix_window = 2e6,
              ens_db = ensDb_v100)
#loc2 <- link_LD(loc2, token = "9c40685f2cfb")

candQtlSub1 <- qtl1[qtl1$pheID == cand, ]
candQtlSub1 <- merge(candQtlSub1, ld2[, c(6,7)], by.x = "rsid", by.y = "SNP_B")
candQtlSub1$logP <- -log10(candQtlSub1$p)
candQtlSub1$ld <- candQtlSub1$R2
loc1 <- locus(data = candQtlSub1, gene = candName, fix_window = 2e6,
              ens_db = ensDb_v100)
#loc1 <- link_LD(loc1, token = "9c40685f2cfb")



locus_plot(loc1, labels = c("index"))
locus_plot(loc2, labels = c("index")) 


ylim_combined <- range(c(-log10(loc1$data$p), -log10(loc2$data$p)))


#png(paste0("resPlots/Fig2/", candName, "_locusplot2.png"), width = 3, height = 3.5, units = "in", res = 300)
pdf(paste0("resPlots/Fig2/", candName, "_locusplot2.pdf"), width = 2.8, height = 3)
# set up layered plot with 2 plots & a gene track; store old par() settings
oldpar <- set_layers(2)
scatter_plot(loc1, xticks = FALSE, labels = c(), ylim = ylim_combined)
scatter_plot(loc2, xticks = FALSE, labels = c("index"), ylim = ylim_combined)
genetracks(loc1, highlight = candName, blanks = c("hide"), filter_gene_name = c(candName, "GARNL3", "NIBAN2", "DPM2", "SET", "TTC16",
                                                                                "RPL12", "CRAT"),
           filter_gene_biotype = 'protein_coding', maxrows = 2)
par(oldpar)  # revert par() settings
dev.off()




