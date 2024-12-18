# This script takes three arguments
# - Ensemble gene ID
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
  make_option(c("-i", "--candGeneID"), type="character", default=NULL, 
              help="Ensemble gene ID. eg ENSG00000050820", metavar="character"),
  make_option(c("-c", "--coordID"), type="character", default=NULL, 
              help="coordID of SNP. eg 9:61667843:G:C", metavar="character"),
  make_option(c("-r", "--rsID"), type="character", default="toAdd", 
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

if (is.null(opt$coordID)){
  print_help(opt_parser)
  stop("Please provide coord ID.", call.=FALSE)
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
cd <- paste0(opt$coordID)
#cd <- "9:61667843:G:C"
rs <- paste0(opt$rsID)
rs <- ifelse(rs == "toAdd" , cd, rs)
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
  geom_jitter(colour="darkgrey", position=position_jitter(width=0.2), size = 0.5) +
  geom_boxplot(outlier.size=0, alpha=0.5, linewidth = 0.05, fill="blue") + 
  labs(x = rs, y = candName, size=7)+
  scale_x_discrete(labels=c(hom1,het,hom2))+
  theme_minimal() +
theme(axis.text.x = element_text(size=7),
      axis.text.y = element_text(size=7),
      axis.title.x = element_text(size=7),
      axis.title.y = element_text(size=7, face="bold"))



ggsave(paste0(pp, "_", cand, "_vs_", cd, "_Plot.pdf" ), dfPlot,
       scale = 1,
       path = paste0("resPlots/Fig2"),
       width = 40,
       height = 35,
       units = "mm",
       dpi = 300,
       limitsize = TRUE)
