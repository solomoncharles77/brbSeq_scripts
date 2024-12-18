library(tidyverse)
library(data.table)
library(biomaRt)


# Functions ---------------------------------------------------------------

getBiotypeAll <- function(biotypeDF){
  biotypeDF <- data.frame(table(biotypeDF$gene_biotype))
  biotypeDF <- biotypeDF[order(biotypeDF$Freq, decreasing = T), ]
  colnames(biotypeDF) <- c("Gene_Biotype", "Frequency")
  return(biotypeDF)
}

getBiotypeSummary <- function(biotypeDF){
  biotypeDF %>% 
    filter(str_detect(Gene_Biotype, "pseudogene")) -> pseudoDF
  biotypeDF %>% 
    filter(!str_detect(Gene_Biotype, "pseudogene|protein_coding|lncRNA")) -> othersDF
  
  biotypeDF %>% 
    add_row(Gene_Biotype = "toAdd", Frequency = sum(pseudoDF$Frequency)) %>% 
    add_row(Gene_Biotype = "Others", Frequency = sum(othersDF$Frequency)) %>% 
    filter(str_detect(Gene_Biotype, "toAdd|protein_coding|lncRNA|Others")) %>% 
    add_row(Gene_Biotype = "Total", Frequency = sum(biotypeDF$Frequency)) %>% 
    mutate(Gene_Biotype = str_replace(Gene_Biotype, "toAdd", "Pseudogenes")) -> biotypeDF1
  
  return(biotypeDF1)
}

qtlHeader <- read.delim("brbSeq_scripts/QTLtools_permute_header.txt", header = F)
colIDs <- qtlHeader$V2
colIDs <- sub(" \\| ve_by_pc1 \\| n_phe_in_grp", "", colIDs)
colIDs <- sub("phe_id \\| ", "", colIDs)

# Get biotype for all tested genes quantified -----------------------
Splice <- data.frame(fread("testSeq/readCount/rawReadCount.csv"))
ensembl <- useEnsembl(biomart = "ensembl",
                      dataset = "hsapiens_gene_ensembl")
ensIDsAnnot <-getBM(
  attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'),
  filters=c('ensembl_gene_id'),
  values = Splice$V1,
  mart = ensembl, uniqueRows = TRUE)

SpliceAnot <- merge(ensIDsAnnot, Splice, by.x = "ensembl_gene_id", by.y = "V1", all.y = TRUE)
SpliceAnot[c("gene_biotype")][is.na(SpliceAnot[c("gene_biotype")])] <- "NA"
SpliceAnot <- SpliceAnot[, 1:3]

# Add TNFa input Sqtl biotype --------------------------------------------------------
tnfa <- data.frame(fread("phenoFiles/LeafcutterSplicing_TNFA_qtlTools_Ready_IntronAnnotation.txt"))
tnfa <- tnfa[!duplicated(tnfa$geneID), ]
tnfaSpliceAnot <- merge(SpliceAnot, tnfa, by.x = "ensembl_gene_id", by.y = "geneID")
tnfaSpliceAnot <- tnfaSpliceAnot[, 1:3]
tnfaSpliceBiotype_all <- getBiotypeAll(tnfaSpliceAnot)
tnfaSpliceBiotype_Summ <- getBiotypeSummary(tnfaSpliceBiotype_all)
# write.csv(tnfaSpliceBiotype_all, "resFiles/tnfaBiotype_of_Sqtl_tested_genes_All.csv", row.names = F)
# write.csv(tnfaSpliceBiotype_Summ, "resFiles/tnfaBiotype_of_Sqtl_tested_genes_Summary.csv", row.names = F)

# Add TNFa Sqtl eGene biotype --------------------------------------------------------
tnfaSqtl <- data.frame(fread("leafcutterSplicing_TNFA_cisEQTL/leafcutterSplicing_TNFA_cisEQTL_permute_ALL.txt.gz"))
colnames(tnfaSqtl) <- paste0("tnfa_", colIDs)
tnfaSqtl <- tnfaSqtl[!duplicated(tnfaSqtl$tnfa_grp_id), ]
tnfaSqtl <- tnfaSqtl[tnfaSqtl$tnfa_adj_emp_pval < 0.05, ]
tnfaSqtl <- tnfaSqtl[order(tnfaSqtl$tnfa_adj_emp_pval), ]
tnfaSqtl <- tnfaSqtl[!duplicated(tnfaSqtl$tnfa_phe_id), ]

tnfaSqtlAnot <- merge(SpliceAnot, tnfaSqtl, by.x = "ensembl_gene_id", by.y = "tnfa_grp_id")
tnfaSqtlAnot <- tnfaSqtlAnot[, 1:3]
tnfaSqtlBiotype_all <- getBiotypeAll(tnfaSqtlAnot)
tnfaSqtlBiotype_Summ <- getBiotypeSummary(tnfaSqtlBiotype_all)


# Add UNTR input Sqtl biotype --------------------------------------------------------
untr <- data.frame(fread("phenoFiles/LeafcutterSplicing_UNTR_qtlTools_Ready_IntronAnnotation.txt"))
untr <- untr[!duplicated(untr$geneID), ]
untrSpliceAnot <- merge(SpliceAnot, untr, by.x = "ensembl_gene_id", by.y = "geneID")
untrSpliceAnot <- untrSpliceAnot[, 1:3]
untrSpliceBiotype_all <- getBiotypeAll(untrSpliceAnot)
untrSpliceBiotype_Summ <- getBiotypeSummary(untrSpliceBiotype_all)
#write.csv(untrSpliceBiotype_all, "resFiles/untrBiotype_of_Sqtl_tested_genes_All.csv", row.names = F)
#write.csv(untrSpliceBiotype_Summ, "resFiles/untrBiotype_of_Sqtl_tested_genes_Summary.csv", row.names = F)

# Add TNFa Sqtl eGene biotype --------------------------------------------------------
untrSqtl <- data.frame(fread("leafcutterSplicing_UNTR_cisEQTL/leafcutterSplicing_UNTR_cisEQTL_permute_ALL.txt.gz"))
colnames(untrSqtl) <- paste0("untr_", colIDs)
untrSqtl <- untrSqtl[!duplicated(untrSqtl$untr_grp_id), ]
untrSqtl <- untrSqtl[untrSqtl$untr_adj_emp_pval < 0.05, ]
untrSqtl <- untrSqtl[order(untrSqtl$untr_adj_emp_pval), ]
untrSqtl <- untrSqtl[!duplicated(untrSqtl$untr_phe_id), ]

untrSqtlAnot <- merge(SpliceAnot, untrSqtl, by.x = "ensembl_gene_id", by.y = "untr_grp_id")
untrSqtlAnot <- untrSqtlAnot[, 1:3]
untrSqtlBiotype_all <- getBiotypeAll(untrSqtlAnot)
untrSqtlBiotype_Summ <- getBiotypeSummary(untrSqtlBiotype_all)






# Plot Naive and TNFa input gene Spliceession biotype barplot -------------------------------------------

untrTnfaSplice <- merge(untrSpliceBiotype_Summ, tnfaSpliceBiotype_Summ, by = "Gene_Biotype")
colnames(untrTnfaSplice) <- c("Gene_Biotype", "Naive", "TNFa")

# Convert data from wide to long format for ggplot
untrTnfaSplice <- untrTnfaSplice %>%
  pivot_longer(cols = c(Naive, TNFa), 
               names_to = "Frequency_Type", 
               values_to = "Frequency") %>%
  mutate(Gene_Biotype = fct_reorder(Gene_Biotype, Frequency))

# Create bar plot using ggplot2
ggplot(untrTnfaSplice, aes(x = Gene_Biotype, y = Frequency, fill = Frequency_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(title = "Spliced Gene Biotypes",
       x = "Gene Biotype",
       y = "Frequency",
       fill = "Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  -> untrTnfaSplicePlot
untrTnfaSplicePlot

ggsave("naiveTnfaSplicePlot.png", plot = untrTnfaSplicePlot,
       path = "resPlots", height = 4.73, width = 7.45, units = "in", bg = "white")

# Plot Naive and TNFa Sqtl eGene biotype barplot -------------------------------------------

untrTnfaSqtl <- merge(untrSqtlBiotype_Summ, tnfaSqtlBiotype_Summ, by = "Gene_Biotype")
colnames(untrTnfaSqtl) <- c("Gene_Biotype", "Naive", "TNFa")

# Convert data from wide to long format for ggplot
untrTnfaSqtl <- untrTnfaSqtl %>%
  pivot_longer(cols = c(Naive, TNFa), 
               names_to = "Frequency_Type", 
               values_to = "Frequency") %>%
  mutate(Gene_Biotype = fct_reorder(Gene_Biotype, Frequency))

# Create bar plot using ggplot2
ggplot(untrTnfaSqtl, aes(x = Gene_Biotype, y = Frequency, fill = Frequency_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(title = "sGene Biotypes",
       x = "Gene Biotype",
       y = "Frequency",
       fill = "Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))   -> untrTnfaSqtlPlot
untrTnfaSqtlPlot

ggsave("naiveTnfaSqtlPlot.png", plot = untrTnfaSqtlPlot,
       path = "resPlots", height = 4.73, width = 7.45, units = "in", bg = "white")
