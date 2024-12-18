# This script takes three arguments
# - Ensemble gene ID
# - rsID of SNP
# - SNP Coord
# This script outputs a hg38 version of the gwas summary statistics file 


# Set global options and load libraries -----------------------------------
options(scipen = 999, "warnPartialMatchDollar"=TRUE)
# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))


# List options
option_list = list(
  make_option(c("-c", "--candGeneID"), type="character", default=NULL, 
              help="Ensemble gene ID. eg ENSG00000050820", metavar="character"),
  make_option(c("-r", "--rsID"), type="character", default=NULL, 
              help="Rs ID of SNP. eg rs11648176", metavar="character"),
  make_option(c("-e", "--expr"), type="character", default=NULL, 
              help="Path to gene expression file", metavar="character"),
  make_option(c("-g", "--geno"), type="character", default=NULL, 
              help="Path to genotype file in dosage format. Just the file prefix e.g. genoFiles/brbSeqGeno_QCed", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default = NULL, 
              help="Prefix for output file.", metavar="character")
); 

opt_parser <-  OptionParser(option_list=option_list);
opt <-  parse_args(opt_parser);

# Data input error messages

if (is.null(opt$candGeneID)){
  print_help(opt_parser)
  stop("Please provide gene ID.", call.=FALSE)
}

if (is.null(opt$rsID)){
  print_help(opt_parser)
  stop("Please provide rs ID.", call.=FALSE)
}

if (is.null(opt$expr)){
  print_help(opt_parser)
  stop("Please provide path to normalized gene expression matrix.", call.=FALSE)
}

if (is.null(opt$expr)){
  print_help(opt_parser)
  stop("Please provide path to genotype file in dosage format.", call.=FALSE)
}

if (is.null(opt$prefix)){
  print_help(opt_parser)
  stop("Please provide prefix that identifies the pheno data", call.=FALSE)
}

# Assign input variables --------------------------------------------------

# Import gene expression --------------------------------------------------
exprFilename <- paste0(opt$expr)
#exprFilename <- "phenoFiles/normGeneExpr_TNFA_qtlTools_Ready.bed.gz"
genoFilename <- paste0(opt$geno)
#genoFilename <- "genoFiles/brbSeqGeno_QCed"
cand <- paste0(opt$candGeneID)
#cand <- "ENSG00000182511"
rs <- paste0(opt$rsID)
#rs <- "rs2521501"
pp <- paste0(opt$prefix)
#pp <- "geneExprTNFA"


# Import data -------------------------------------------------------------
expr <- data.frame(fread(exprFilename))
expr <- expr[, c(5,7:ncol(expr))]
annot <- data.frame(fread("exprFiles/geneAnnot.csv"))

# Extract expr for cand gene
candexpr <- expr[expr$gene %in% cand, ]
candName <- annot[annot$ensembl_gene_id %in% cand, 2]

# Fetch rsID details
rsFile <- paste0("get_rsID_Coord.txt")
write.table(rs, file = rsFile, quote = F, row.names = F, col.names = F)
system(paste0("grep -wFf ", rsFile, " /scratch/vasccell/cs806/colocalization/dbSNP/hg38.snp151_All.bed | cut -f4,3,10 > rsLookup.txt"))
rsCoord <- read.table("rsLookup.txt")
colnames(rsCoord) <- c("bp", "rsID", "coord")

# TODO
# I need some if commands around here to handle multiple coord cases.
#pos <- as.numeric(rsCoord$bp)

# if (nrow(rsCoord) > 1) {
#   cat("\n", rs, "returns more than one coordinate. \n")
#   cat("\n Checking which coordinate is present in genotype data \n")
#   for (rr in nrow(rsCoord)) {
#     nn = rsCoord$coord[rr]
#     nnLoc <- system(paste0("grep -n ", nn, " ", genoFilename, "_Geno_snpsloc.txt"), intern = T)
#     if(nnLoc == NULL)
#   }
#   cd <- rsCoord$coord
#   return(cd)
#   
# }



cd <- rsCoord$coord[1]
#cd <- rsCoord$coord[2]

a1 <- unlist(strsplit(cd, split = ":"))[3]
a2 <- unlist(strsplit(cd, split = ":"))[4]
hom1 <- paste0(a1, a1)
het <- paste0(a1, a2)
hom2 <- paste0(a2, a2)


# Extract geno for cand SNP
exprLoc <- system(paste0("grep -n ", cd, " ", genoFilename, "_Geno_snpsloc.txt"), intern = T)
toSkip <- as.numeric(gsub(":.*", "", exprLoc))
toSkip <- toSkip - 1
geno <- data.frame(fread(paste0(genoFilename, "_Geno.txt"), skip = toSkip, nrows = 1, col.names=scan(paste0(genoFilename, "_MatrixHeader.txt"), what ="", sep = "\t", quiet = TRUE)))


# Harmonize geno and expr
colnames(geno) <- gsub(".*_S", "S", colnames(geno))
genoSamps <- colnames(geno)[c(colnames(geno) %in% colnames(candexpr))]

candexpr <- candexpr[, genoSamps]
candexpr <- data.frame(Gene = cand, candexpr)

df <- t(rbind(geno[-1], candexpr[-1]))
colnames(df) = c("snp", "gene_expr")
df <- data.frame(df)
df <- df[complete.cases(df), ]
df$snp = as.factor(df$snp)

dflm = lm(df[,"gene_expr"] ~ as.numeric(df[,"snp"]))

dfPlot <- ggplot(df, aes(snp, gene_expr)) +
  geom_jitter(colour="darkgrey", position=position_jitter(width=0.25)) +
  geom_boxplot(outlier.size=0, alpha=0.6, fill="lightgreen") + 
  labs(x = rs, y = candName)+
  scale_x_discrete(labels=c(hom1,het,hom2))+
  theme_bw()
# theme(axis.text.x = element_text(face="bold", color="#993333", 
#                                  size=14, angle=45))

ggsave(paste0(pp, "_", cand, "_vs_", rs, "_Plot.png" ), dfPlot,
       scale = 1,
       path = paste0("resPlots/exprGenoPlot"),
       width = 80,
       height = 70,
       units = "mm",
       dpi = 300,
       limitsize = TRUE)
