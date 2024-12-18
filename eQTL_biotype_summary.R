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
expr <- data.frame(fread("testSeq/readCount/rawReadCount.csv"))
expr[1:10, 1:10]

ensembl <- useEnsembl(biomart = "ensembl",
                      dataset = "hsapiens_gene_ensembl")
ensIDsAnnot <-getBM(
  attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'),
  filters=c('ensembl_gene_id'),
  values = expr$V1,
  mart = ensembl, uniqueRows = TRUE)

exprAnot <- merge(ensIDsAnnot, expr, by.x = "ensembl_gene_id", by.y = "V1", all.y = TRUE)
exprAnot[c("gene_biotype")][is.na(exprAnot[c("gene_biotype")])] <- "NA"
exprAnot <- exprAnot[, 1:3]

# Add TNFa input eQTL biotype --------------------------------------------------------
tnfa <- data.frame(fread("phenoFiles/normGeneExpr_TNFA_qtlTools_Ready.bed.gz"))
tnfaExprAnot <- merge(exprAnot, tnfa, by.x = "ensembl_gene_id", by.y = "gene")
tnfaExprAnot <- tnfaExprAnot[, 1:3]
tnfaExprBiotype_all <- getBiotypeAll(tnfaExprAnot)
tnfaExprBiotype_Summ <- getBiotypeSummary(tnfaExprBiotype_all)
write.csv(tnfaExprBiotype_all, "resFiles/tnfaBiotype_of_eQTL_tested_genes_All.csv", row.names = F)
write.csv(tnfaExprBiotype_Summ, "resFiles/tnfaBiotype_of_eQTL_tested_genes_Summary.csv", row.names = F)

# Add TNFa eQTL eGene biotype --------------------------------------------------------
tnfaEqtl <- data.frame(fread("geneExprTNFA_brbSeq_cisEQTL/geneExprTNFA_brbSeq_cisEQTL_permute_ALL.txt.gz"))
colnames(tnfaEqtl) <- paste0("tnfa_", colIDs)
tnfaEqtl <- tnfaEqtl[tnfaEqtl$tnfa_adj_emp_pval < 0.05, ]
tnfaEqtl <- tnfaEqtl[order(tnfaEqtl$tnfa_adj_emp_pval), ]
tnfaEqtl <- tnfaEqtl[!duplicated(tnfaEqtl$tnfa_phe_id), ]

tnfaEqtlAnot <- merge(exprAnot, tnfaEqtl, by.x = "ensembl_gene_id", by.y = "tnfa_phe_id")
tnfaEqtlAnot <- tnfaEqtlAnot[, 1:3]
tnfaEqtlBiotype_all <- getBiotypeAll(tnfaEqtlAnot)
tnfaEqtlBiotype_Summ <- getBiotypeSummary(tnfaEqtlBiotype_all)


# Add UNTR input eQTL biotype --------------------------------------------------------
untr <- data.frame(fread("phenoFiles/normGeneExpr_UNTR_qtlTools_Ready.bed.gz"))
untrExprAnot <- merge(exprAnot, untr, by.x = "ensembl_gene_id", by.y = "gene")
untrExprAnot <- untrExprAnot[, 1:3]
untrExprBiotype_all <- getBiotypeAll(untrExprAnot)
untrExprBiotype_Summ <- getBiotypeSummary(untrExprBiotype_all)
write.csv(untrExprBiotype_all, "resFiles/untrBiotype_of_eQTL_tested_genes_All.csv", row.names = F)
write.csv(untrExprBiotype_Summ, "resFiles/untrBiotype_of_eQTL_tested_genes_Summary.csv", row.names = F)

# Add TNFa eQTL eGene biotype --------------------------------------------------------
qtlHeader <- read.delim("brbSeq_scripts/QTLtools_permute_header.txt", header = F)
colIDs <- qtlHeader$V2
colIDs <- sub(" \\| ve_by_pc1 \\| n_phe_in_grp", "", colIDs)
colIDs <- sub("phe_id \\| ", "", colIDs)

untrEqtl <- data.frame(fread("geneExprUNTR_brbSeq_cisEQTL/geneExprUNTR_brbSeq_cisEQTL_permute_ALL.txt.gz"))
colnames(untrEqtl) <- paste0("untr_", colIDs)
untrEqtl <- untrEqtl[untrEqtl$untr_adj_emp_pval < 0.05, ]
untrEqtl <- untrEqtl[order(untrEqtl$untr_adj_emp_pval), ]
untrEqtl <- untrEqtl[!duplicated(untrEqtl$untr_phe_id), ]

untrEqtlAnot <- merge(exprAnot, untrEqtl, by.x = "ensembl_gene_id", by.y = "untr_phe_id")
untrEqtlAnot <- untrEqtlAnot[, 1:3]
untrEqtlBiotype_all <- getBiotypeAll(untrEqtlAnot)
untrEqtlBiotype_Summ <- getBiotypeSummary(untrEqtlBiotype_all)






# Plot Naive and TNFa input gene expression biotype barplot -------------------------------------------

untrTnfaExpr <- merge(untrExprBiotype_Summ, tnfaExprBiotype_Summ, by = "Gene_Biotype")
colnames(untrTnfaExpr) <- c("Gene_Biotype", "Naive", "TNFa")

# Convert data from wide to long format for ggplot
untrTnfaExpr <- untrTnfaExpr %>%
  pivot_longer(cols = c(Naive, TNFa), 
               names_to = "Frequency_Type", 
               values_to = "Frequency") %>%
  mutate(Gene_Biotype = fct_reorder(Gene_Biotype, Frequency))

# Create bar plot using ggplot2
ggplot(untrTnfaExpr, aes(x = Gene_Biotype, y = Frequency, fill = Frequency_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(title = "Gene Expression Biotypes",
       x = "Gene Biotype",
       y = "Frequency",
       fill = "Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> untrTnfaExprPlot
untrTnfaExprPlot

ggsave("naiveTnfaExprPlot.png", plot = untrTnfaExprPlot,
       path = "resPlots", height = 4.73, width = 7.45, units = "in", bg = "white")

# Plot Naive and TNFa eQTL eGene biotype barplot -------------------------------------------

untrTnfaEqtl <- merge(untrEqtlBiotype_Summ, tnfaEqtlBiotype_Summ, by = "Gene_Biotype")
colnames(untrTnfaEqtl) <- c("Gene_Biotype", "Naive", "TNFa")

# Convert data from wide to long format for ggplot
untrTnfaEqtl <- untrTnfaEqtl %>%
  pivot_longer(cols = c(Naive, TNFa), 
               names_to = "Frequency_Type", 
               values_to = "Frequency") %>%
  mutate(Gene_Biotype = fct_reorder(Gene_Biotype, Frequency))

# Create bar plot using ggplot2
ggplot(untrTnfaEqtl, aes(x = Gene_Biotype, y = Frequency, fill = Frequency_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(title = "eGene Biotypes",
       x = "Gene Biotype",
       y = "Frequency",
       fill = "Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  -> untrTnfaEqtlPlot
untrTnfaEqtlPlot

ggsave("naiveTnfaEqtlPlot.png", plot = untrTnfaEqtlPlot,
       path = "resPlots", height = 4.73, width = 7.45, units = "in", bg = "white")
