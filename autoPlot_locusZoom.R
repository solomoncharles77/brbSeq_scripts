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
suppressPackageStartupMessages(library(AnnotationHub))
suppressPackageStartupMessages(library(locuszoomr))


# List options
option_list = list(
  make_option(c("-c", "--candsFile"), type="character", default=NULL, 
              help="Ensembl gene ID. eg ENSG00000050820", metavar="character"),
  make_option(c("-n", "--geneName"), type="character", default=NULL, 
              help="External gene name eg FES", metavar="character")
); 

opt_parser <-  OptionParser(option_list=option_list);
opt <-  parse_args(opt_parser);

# Data input error messages

if (is.null(opt$candsFile)){
  print_help(opt_parser)
  stop("Please provide gene ID.", call.=FALSE)
}


# Import data -------------------------------------------------------------
ah <- AnnotationHub()
query(ah, c("EnsDb", "v100", "Homo sapiens"))
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


# Assign variables --------------------------------------------------------
candsPath <- paste0(opt$candsFile)
#candsPath <- "resFiles/trsCands.txt"
cands <- read.delim(candsPath, header = F)


#candName <- if (is.null(candName)) cand else candName

# for (x in 1:nrow(cands)) {
#   
#   cand <- cands$V1[x]
#   candName <- cands$V2[x]
#   
#   if (is.null(candName) || is.na(candName) || trimws(candName) == "") {
#     candName <- cand
#   }
#   
#   # Format data -------------------------------------------------------------
#   candQtlSub1 <- qtl1[qtl1$pheID == cand, ]
#   loc1 <- locus(data = candQtlSub1, gene = candName, fix_window = 2e6,
#                 ens_db = ensDb_v100)
#   #loc1 <- link_LD(loc1, token = "9c40685f2cfb")
#   
#   candQtlSub2 <- qtl2[qtl2$pheID == cand, ]
#   loc2 <- locus(data = candQtlSub2, gene = candName, fix_window = 2e6,
#                 ens_db = ensDb_v100)
#   #loc2 <- link_LD(loc2, token = "9c40685f2cfb")
#   
#   
#   ylim_combined <- range(c(-log10(loc1$data$p), -log10(loc2$data$p)))
#   
#   #pdf("resPlots/myplot2.pdf", width = 6, height = 8)
#   png(paste0("resPlots/locusPlots/", candName, "_locusplot.png"), width = 6, height = 8, units = "in", res = 300)
#   # set up layered plot with 2 plots & a gene track; store old par() settings
#   oldpar <- set_layers(2)
#   scatter_plot(loc1, xticks = FALSE, labels = c("index"), ylim = ylim_combined)
#   scatter_plot(loc2, xticks = FALSE, labels = c("index"), ylim = ylim_combined)
#   genetracks(loc1, highlight = candName, blanks = c("hide"))
#   par(oldpar)  # revert par() settings
#   dev.off()
#   
# }

# Loop through the candidates
for (x in 1:nrow(cands)) {
  
  cand <- cands$V1[x]
  candName <- cands$V2[x]
  
  if (is.null(candName) || is.na(candName) || trimws(candName) == "") {
    candName <- cand
  }
  
  # Wrap the loop body in a try block
  result <- try({
    # Format data -------------------------------------------------------------
    candQtlSub1 <- qtl1[qtl1$pheID == cand, ]
    loc1 <- locus(data = candQtlSub1, gene = candName, fix_window = 2e6,
                  ens_db = ensDb_v100)
    
    candQtlSub2 <- qtl2[qtl2$pheID == cand, ]
    loc2 <- locus(data = candQtlSub2, gene = candName, fix_window = 2e6,
                  ens_db = ensDb_v100)
    
    ylim_combined <- range(c(-log10(loc1$data$p), -log10(loc2$data$p)))
    
    # Save plot
    png(paste0("resPlots/locusPlots/", candName, "_locusplot.png"), width = 6, height = 8, units = "in", res = 300)
    # set up layered plot with 2 plots & a gene track; store old par() settings
    oldpar <- set_layers(2)
    scatter_plot(loc1, xticks = FALSE, labels = c("index"), ylim = ylim_combined)
    scatter_plot(loc2, xticks = FALSE, labels = c("index"), ylim = ylim_combined)
    genetracks(loc1, highlight = candName, blanks = c("hide"))
    par(oldpar)  # revert par() settings
    dev.off()
    
  }, silent = TRUE)  # This makes sure the error message is not printed automatically
  
  # Check if an error occurred
  if (inherits(result, "try-error")) {
    message("Error encountered for candidate: ", cand, " with gene: ", candName)
    next  # Skip to the next iteration
  }
  
  # If no error, print the result (or any other desired actions)
  # print(result)
}

