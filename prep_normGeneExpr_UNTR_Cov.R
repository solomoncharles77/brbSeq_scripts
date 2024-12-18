
sex <- read.table("genoFiles/BRB-seq_Pilot_Study_All_SNPs_No_Filtering.fam")
sex <- data.frame(t(sex[, c(2,5)]))
colnames(sex) <- sex[1, ]
sex <- sex[-1, ]
sex <- sex[, order(colnames(sex))]
sex <- cbind(SampleID = "Sex", sex)

genoPCA <- read.table(paste0("covFiles/brbSeqGeno_genoPCs.txt"), check.names=FALSE)
genoPCA <- genoPCA[1:3, ]
phenoPCA <- read.table(paste0("covFiles/normGeneExpr_UNTR_qtlTools_Ready.pca"), header = T)
phenoPCA <- phenoPCA[1:40, ]



cov <- rbind(sex, genoPCA, phenoPCA)

write.table(cov, paste0("covFiles/normGeneExpr_UNTR_geno3pc_pheno40pc.txt"), sep = "\t", row.names = F, quote = F)

cat("Covariates file Ready  \n")
