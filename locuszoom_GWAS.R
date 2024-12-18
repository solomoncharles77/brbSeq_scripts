
library(data.table)
library(locuszoomr)
library(AnnotationHub)

ah <- AnnotationHub()
query(ah, c("EnsDb", "v100", "Homo sapiens"))
#ensDb_v106 <- ah[["AH100643"]]
ensDb_v100 <- ah[["AH79689"]]


cand <- "ENSG00000106367"
candName <- "AP1S1"

geneCoords <- data.frame(fread("../exprPhenoData/VSMC_Gene_Coords.txt"))
candCoord <- geneCoords[geneCoords$geneID %in% cand, ]
candCoord$up1mb <- candCoord$Start + 1000000
candCoord$down1mb <- candCoord$Start - 1000000

gwas1 <- data.frame(fread("../colocalization/cleanGWAS_Summary_Stats/GWAS_Coronary-Artery-Disease_Aragam_2022_NatGenet_hg38.txt"))

gwas1 <- gwas1[, c("permID", "chr", "hg38_bp", "REF", "ALT", "pvalue", "beta", "se")]
colnames(gwas1) <- c("rsid", "chrom", "pos", "other_allele", "effect_allele", "p", "beta", "se")

gwas1Sub <- gwas1[gwas1$chrom == candCoord$Chr & gwas1$pos >= candCoord$down1mb & gwas1$pos < candCoord$up1mb, ]


gwas1Sub <- gwas1Sub[order(gwas1Sub$p, gwas1Sub$rsid), ]
topSNP <- gwas1Sub[1, 1]

snpListFile <- paste0(candName, "_SNPList.txt")
topSNPFile <- paste0(candName, "_topSNP.txt")
snpLDFile <- paste0(candName, "_SNPList_LD.txt")
write.table(data.frame(gwas1Sub$rsid), file = snpListFile,  row.names = F, col.names = F, quote = F)
write.table(data.frame(topSNP), file = topSNPFile,  row.names = F, col.names = F, quote = F)
system(paste0("plink --bfile /scratch/vasccell/cs806/colocalization/1000Genome/g1000_eur --extract ", snpListFile, " --ld-window-r2 0.00000001 --r2  --ld-window-kb 2000000 --ld-window 1000000 --ld-snp-list ", topSNPFile, " --threads 6 --out ", snpLDFile))
gwasLD <- read.table(paste0(snpLDFile, ".ld"), header = T)

gwas1Sub <- merge(gwas1Sub, gwasLD[, c(6,7)], by.x = "rsid", by.y = "SNP_B")
gwas1Sub$logP <- -log10(gwas1Sub$p)
gwas1Sub$ld <- gwas1Sub$R2
gwas1Loc <- locus(data = gwas1Sub, gene = candName, fix_window = 2e6,
              ens_db = ensDb_v100)

locus_plot(gwas1Loc, labels = c("index"))




qtl1 <- data.frame(fread("geneExprUNTR_brbSeq_cisEQTL/geneExprUNTR_brbSeq_cisEQTL_nominal_ALL_pval0.05_rsID.txt.gz",
                         select = c(7,25,11,12,19,22,23,14,16,17,15),
                         col.names = c("pheID", "GeneName", "chrom", "pos", "rsid",
                                       "other_allele", "effect_allele", "p", "beta", "se", "r2")))
candQtlSub1 <- qtl1[qtl1$pheID == cand, ]
candQtlSub1 <- candQtlSub1[order(candQtlSub1$p, candQtlSub1$rsid), ]
topSNP <- candQtlSub1[1, 5]

snpListFile <- paste0(candName, "_SNPList.txt")
topSNPFile <- paste0(candName, "_topSNP.txt")
snpLDFile <- paste0(candName, "_SNPList_LD.txt")
write.table(data.frame(candQtlSub1$rsid), file = snpListFile,  row.names = F, col.names = F, quote = F)
write.table(data.frame(topSNP), file = topSNPFile,  row.names = F, col.names = F, quote = F)
system(paste0("plink --bfile /scratch/vasccell/cs806/colocalization/1000Genome/g1000_eur --extract ", snpListFile, " --ld-window-r2 0.00000001 --r2  --ld-window-kb 2000000 --ld-window 1000000 --ld-snp-list ", topSNPFile, " --threads 6 --out ", snpLDFile))
qtl1LD <- read.table(paste0(snpLDFile, ".ld"), header = T)

candQtlSub1 <- merge(candQtlSub1, qtl1LD[, c(6,7)], by.x = "rsid", by.y = "SNP_B")
candQtlSub1$logP <- -log10(candQtlSub1$p)
candQtlSub1$ld <- candQtlSub1$R2
qtl1Loc <- locus(data = candQtlSub1, gene = candName, fix_window = 2e6,
              ens_db = ensDb_v100)
#loc1 <- link_LD(loc1, token = "9c40685f2cfb")
locus_plot(qtl1Loc, labels = c("index"))




qtl2 <- data.frame(fread("geneExprTNFA_brbSeq_cisEQTL/geneExprTNFA_brbSeq_cisEQTL_nominal_ALL_pval0.05_rsID.txt.gz",
                         select = c(7,25,11,12,19,22,23,14,16,17,15),
                         col.names = c("pheID", "GeneName", "chrom", "pos", "rsid",
                                       "other_allele", "effect_allele", "p", "beta", "se", "r2")))
gc()
candQtlSub2 <- qtl2[qtl2$pheID == cand, ]
candQtlSub2 <- merge(candQtlSub2, ld[, c(6,7)], by.x = "rsid", by.y = "SNP_B")
candQtlSub2$logP <- -log10(candQtlSub2$p)
candQtlSub2$ld <- candQtlSub2$R2
loc2 <- locus(data = candQtlSub2, gene = candName, fix_window = 2e6,
              ens_db = ensDb_v100)
#loc2 <- link_LD(loc2, token = "9c40685f2cfb")
locus_plot(loc2, labels = c("index")) 


ylim_combined <- range(c(-log10(gwas1Loc$data$p), -log10(qtl1Loc$data$p)))

#png(paste0("resPlots/Fig2/", candName, "_locusplot2.png"), width = 3, height = 3.5, units = "in", res = 300)
pdf(paste0("resPlots/Fig2/", candName, "_locusplot_test.pdf"), width = 2.8, height = 3)
# set up layered plot with 2 plots & a gene track; store old par() settings
oldpar <- set_layers(2)
scatter_plot(gwas1Loc, xticks = FALSE, labels = c("index"), ylim = ylim_combined)
scatter_plot(qtl1Loc, xticks = FALSE, labels = c(), ylim = ylim_combined)
genetracks(qtl1Loc, highlight = candName, blanks = c("hide"), filter_gene_name = c(candName, "CUX1", "IFT22", "MUC3A", "AGFG2", "COL26A1",
                                                                                "MYL10"),
           filter_gene_biotype = 'protein_coding', maxrows = 2)
par(oldpar)  # revert par() settings
dev.off()




