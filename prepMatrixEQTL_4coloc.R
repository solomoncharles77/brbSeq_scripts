# This script takes two arguments
# - the path to the eQTL summstats file
# - Path to .frq file for geneotype data used for the eQTL analyses
# This script outputs multiple files used by various colocalization tools
# to the source directory of the summstats file.


# Set global options and load libraries -----------------------------------
options(scipen = 999, "warnPartialMatchDollar"=TRUE)
# Load libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(gtools))


# List options
option_list = list(
  make_option(c("-q", "--qtl"), type="character", default=NULL, 
              help="path to qtl summary statistics", metavar="character"),
  make_option(c("-f", "--freq"), type="character", default=NULL, 
              help="path to plink generated .frq of genotype file", metavar="character")
); 


opt_parser <-  OptionParser(option_list=option_list);
opt <-  parse_args(opt_parser);

# Data input error messages
if (is.null(opt$qtl)){
  print_help(opt_parser)
  stop("Please provide path to summ stats of eqtl.", call.=FALSE)
}

if (is.null(opt$freq)){
  print_help(opt_parser)
  stop("Please provide path to plink generated .frq genotype file.", call.=FALSE)
}


# Assign variables --------------------------------------------------------
qtlFileName <- paste0(opt$qtl)
#qtlFileName <- "geneExpr_cisQTL_3pc_Results/geneExpr_cisQTL_3pc.csv"
frqFileName <- paste0(opt$freq)
#frqFileName <- "huvecImputeGeno/huvecImputeEQTL.frq"

outDir <- dirname(qtlFileName)
# outPrefix <- basename(qtlFileName)
outPrefix <- sub(".csv", "", qtlFileName)
# outPrefix <- paste0(outDir, "/", outPrefix)


# Prep Matrix eQTL output for colocalization ----------------------------------------------------
eQTL <- data.frame(fread(qtlFileName))
frq <- data.frame(fread(frqFileName))
eQTL$se <- eQTL$beta/eQTL$statistic
eQTL$variance <- eQTL$se * eQTL$se
eQTL$snps <- sub("chr", "", eQTL$snps)
frq$position <- sapply(frq$SNP, function(x){unlist(strsplit(x, split = ":"))[2]})
frq$SNP <- sub("chr", "", frq$SNP)
eQTL2 <- merge(eQTL, frq, by.x = "snps", by.y = "SNP")
fwrite(eQTL2, paste0(outPrefix, "_v2.txt.gz"), sep = "\t")
gc()


# Make MAF Beta plot ------------------------------------------------------
eqtl_MBp <- ggplot(data=eQTL2, aes(x=MAF, y=beta)) +
  geom_point()+
  theme_bw()

ggsave(paste0(outPrefix, "_MAF_BETA_Plot.png"), eqtl_MBp,
       path = "mafBetaPlots")

# prep for MR -------------------------------------------------------------
rs <- data.frame(fread("/scratch/vasccell/cs806/genotypeData/rawGeno/Imputed_rsLookup.txt"))
colnames(rs) <- c("position", "rsid", "SNP" )
rsCoord <- merge(frq, rs, by = "SNP")
eQTL3 <- merge(eQTL, rsCoord, by.x = "snps", by.y = "SNP")
colnames(eQTL3)[9:12] <- c("chr", "alt", "ref", "maf")
fwrite(eQTL3, paste0(outPrefix, "_4MR.txt.gz"), sep = "\t")
gc()

# Summarize QTL by chromosome ---------------------------------------------
eQTL_List <- split(eQTL3, eQTL3$chr)
ss <- lapply(eQTL_List, function(df) {
  assoc <- nrow(df)
  rr1 <- range(df$beta)[1]
  rr2 <- range(df$beta)[2]
  mf1 <- range(df$maf)[1]
  mf2 <- range(df$maf)[2]
  kk <- cbind(assoc, rr1, rr2, mf1, mf2)
})

ff <- paste0("chr", names(eQTL_List))

ssDF <- do.call(rbind.data.frame, ss)
ssDF <- cbind(contig = ff, ssDF)
rownames(ssDF) <- ff
ssDF <- ssDF[mixedsort(rownames(ssDF)), ]
ssDF$BETArange <- paste0(round(ssDF$rr1, 2), " to +", round(ssDF$rr2, 2))
ssDF$MAFrange <- paste0(" ", round(ssDF$mf1, 2), " to +", round(ssDF$mf2, 2))
ssDF <- ssDF[, c(1,2,7,8)]

write.csv(ssDF, paste0(outPrefix, "_betaMaf_summary.csv"), row.names = F)

# Prep for SMR HEIDI  ---------------------------------------------------------
#eQTL2 <- data.frame(fread(paste0(outPrefix, "_v2.txt")))
eQTL2 <- eQTL2[, c("snps", "gene", "beta", "statistic", "pvalue", "FDR")]
colnames(eQTL2) <- c("SNP", "gene", "beta", "t-stat", "p-value", "FDR")
eQTL2 <- eQTL2[!duplicated(eQTL2[c("SNP","gene")]), ]
fwrite(eQTL2, paste0(outPrefix, "_mateQTL.txt"), sep = "\t")

system(paste0("/home/c/cs806/SMR/smr_v1.3.1_linux_x86_64_static --eqtl-summary ", paste0(outPrefix, "_mateQTL.txt"),  " --matrix-eqtl-format --make-besd --out ", outDir, "/tempSMR"))

# Restart R to clear memory --
esi <- fread(paste0(outDir, "/tempSMR.esi"))
frq <- frq[, c(1:5,7)]
esi$id  <- 1:nrow(esi)
esi <- merge(esi, frq, by.x = "V2", by.y = "SNP")
esi <- esi[order(esi$id), ]
esi <- esi[, c("CHR", "V2", "V3", "position", "A1", "A2", "MAF")]
fwrite(esi, paste0(outPrefix, ".esi"), sep = "\t", col.names = FALSE)

epi <- fread(paste0(outDir, "/tempSMR.epi"))
genePos <- fread("/scratch/vasccell/cs806/exprPhenoData/VSMC_Gene_Coords_withStrand.txt")
epi$id  <- 1:nrow(epi)
epi <- merge(epi, genePos, by.x = "V2", by.y = "geneID")
epi <- epi[order(epi$id), ]
epi <- epi[, c("Chr", "V2", "V3", "Start", "V2", "strand")]
fwrite(epi, paste0(outPrefix, ".epi"), sep = "\t", col.names = FALSE)

system(paste0("rm ", outDir, "/tempSMR*"))
system(paste0("rm ", outPrefix, "_mateQTL.txt"))


#############################################################################################################################
#############################################################################################################################
