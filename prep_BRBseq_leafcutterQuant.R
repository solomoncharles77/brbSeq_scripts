
# load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library(data.table))

# Import and combine chromosome specific data intron counts for TNFA ---------------
allChr = do.call(rbind, lapply(c(1:22), function(x) 
  read.table(paste0("testSeq/mappedReads/brbSeqJunc/TNFA_juncClus_perind.counts.gz.qqnorm_", x))))

# Get colnames from first line of chromosome 1
con <- file(paste0("testSeq/mappedReads/brbSeqJunc/TNFA_juncClus_perind.counts.gz.qqnorm_1"),"r")
first_line <- readLines(con,n=1)
close(con)
colnames(allChr) <- unlist(strsplit(first_line, split = "\t"))

# Prepare intron matrix, positions, principal components ------------------------
# Prep splicing data
sd <- allChr[, c(5:ncol(allChr))]
names(sd) <-  gsub("_TNFa.Aligned.out.sam.sorted.bam", "", names(sd))
names(sd) <- paste0("S", colnames(sd))
sd <- sd[, order(names(sd))]
sd <- cbind(LCC_ID = allChr$ID, sd)


# Annotate splicing clusters
sp <- allChr[, c(1:4)]
sp <- sp[, c(4, 1:3)]
colnames(sp) <- c("LCC_ID", "Chr", "Start", "End")

sp$Strand <-  sapply(sp$LCC_ID, function(x){unlist(strsplit(x, split = "_"))[3]})
sp$Intron_Cluster_ID <- paste0(sp$Chr, ":", sp$Start, ":", sp$End, ":", sp$Strand)
spBED <- sp[, c("Chr", "Start", "End", "Strand")]
fwrite(spBED, file = paste0("exprFiles/testJuncSpliceIntronCoordinatesBED.txt"), sep = "\t",  col.names = F, row.names = F, quote = F)

bed <- system(paste0("bedtools intersect -a ../exprPhenoData/Homo_sapiens.GRCh38.100.gtf.bed.sorted -b exprFiles/testJuncSpliceIntronCoordinatesBED.txt -wb"), intern = TRUE)
bed <- do.call(rbind, strsplit(bed, '\t'))
bed <- data.frame(bed)
bed$Intron_Cluster_ID <- paste0(bed$X7, ":", bed$X8, ":", bed$X9, ":", bed$X10)
bed <- bed[!duplicated(bed$Intron_Cluster_ID), ]
annot <- merge(bed[, c(4,5,11)], sp,  by = "Intron_Cluster_ID", all.y = TRUE)
annot <- annot[, -c(1)]
colnames(annot)[1:2] <- c("geneID", "geneName")

# Merge annotation with splicing data
annotSD <- merge(annot, sd, by = "LCC_ID")
annotSD <- annotSD[, c(4:6,1:2,7, 8:ncol(annotSD))]
annotSD <- annotSD[order(annotSD$Chr, annotSD$Start), ]
colnames(annotSD)[1:6] <- c("#chr", "start", "end", "feature", "gene", "strand")


# Export files for QTLtools ---------------------------------------------
sFileNamePrefix <- paste0("phenoFiles/LeafcutterSplicing_TNFA_qtlTools_Ready")
sPhenoFileName <- paste0(sFileNamePrefix, ".bed")
sPhenoAnnotFileName <- paste0(sFileNamePrefix, "_IntronAnnotation.txt")

# Export and compress pheno bed file
fwrite(annotSD, file = sPhenoFileName, sep = "\t")
system(paste0("bgzip -f ", sPhenoFileName))
system(paste0("tabix -f ", sPhenoFileName, ".gz"))

# Export annotation
write.table(annot, file = sPhenoAnnotFileName, row.names = FALSE, sep = '\t')
rm(list = ls())

##############################################################################################


# Import and combine chromosome specific data intron counts for UNTR ---------------
allChr = do.call(rbind, lapply(c(1:22), function(x) 
  read.table(paste0("testSeq/mappedReads/brbSeqJunc/UNTR_juncClus_perind.counts.gz.qqnorm_", x))))

# Get colnames from first line of chromosome 1
con <- file(paste0("testSeq/mappedReads/brbSeqJunc/UNTR_juncClus_perind.counts.gz.qqnorm_1"),"r")
first_line <- readLines(con,n=1)
close(con)
colnames(allChr) <- unlist(strsplit(first_line, split = "\t"))

# Prepare intron matrix, positions, principal components ------------------------
# Prep splicing data
sd <- allChr[, c(5:ncol(allChr))]
names(sd) <-  gsub("_untr.Aligned.out.sam.sorted.bam", "", names(sd))
names(sd) <- paste0("S", colnames(sd))
sd <- sd[, order(names(sd))]
sd <- cbind(LCC_ID = allChr$ID, sd)


# Annotate splicing clusters
sp <- allChr[, c(1:4)]
sp <- sp[, c(4, 1:3)]
colnames(sp) <- c("LCC_ID", "Chr", "Start", "End")

sp$Strand <-  sapply(sp$LCC_ID, function(x){unlist(strsplit(x, split = "_"))[3]})
sp$Intron_Cluster_ID <- paste0(sp$Chr, ":", sp$Start, ":", sp$End, ":", sp$Strand)
spBED <- sp[, c("Chr", "Start", "End", "Strand")]
fwrite(spBED, file = paste0("exprFiles/testJuncSpliceIntronCoordinatesBED.txt"), sep = "\t",  col.names = F, row.names = F, quote = F)

bed <- system(paste0("bedtools intersect -a ../exprPhenoData/Homo_sapiens.GRCh38.100.gtf.bed.sorted -b exprFiles/testJuncSpliceIntronCoordinatesBED.txt -wb"), intern = TRUE)
bed <- do.call(rbind, strsplit(bed, '\t'))
bed <- data.frame(bed)
bed$Intron_Cluster_ID <- paste0(bed$X7, ":", bed$X8, ":", bed$X9, ":", bed$X10)
bed <- bed[!duplicated(bed$Intron_Cluster_ID), ]
annot <- merge(bed[, c(4,5,11)], sp,  by = "Intron_Cluster_ID", all.y = TRUE)
annot <- annot[, -c(1)]
colnames(annot)[1:2] <- c("geneID", "geneName")

# Merge annotation with splicing data
annotSD <- merge(annot, sd, by = "LCC_ID")
annotSD <- annotSD[, c(4:6,1:2,7, 8:ncol(annotSD))]
annotSD <- annotSD[order(annotSD$Chr, annotSD$Start), ]
colnames(annotSD)[1:6] <- c("#chr", "start", "end", "feature", "gene", "strand")


# Export files for QTLtools ---------------------------------------------
sFileNamePrefix <- paste0("phenoFiles/LeafcutterSplicing_UNTR_qtlTools_Ready")
sPhenoFileName <- paste0(sFileNamePrefix, ".bed")
sPhenoAnnotFileName <- paste0(sFileNamePrefix, "_IntronAnnotation.txt")

# Export and compress pheno bed file
fwrite(annotSD, file = sPhenoFileName, sep = "\t")
system(paste0("bgzip -f ", sPhenoFileName))
system(paste0("tabix ", sPhenoFileName, ".gz"))

# Export annotation
write.table(annot, file = sPhenoAnnotFileName, row.names = FALSE, sep = '\t')
