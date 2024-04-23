# This script takes one argument -  the prefix name of the plink geno output
aha <- commandArgs(trailingOnly = TRUE)

aha2 <- paste0(aha[1], ".fam")

aha3 <- paste0("covFiles/", aha[2], "_Sex.txt")

# read in plink .fam and extract sex
sex <- read.table(aha2, header = FALSE)
sex <- data.frame(t(sex[, c(2,5)]))
colnames(sex) <- sex[1,]
sex <- sex[-1, ]
sex <- data.frame(Covariates = "Sex", sex)

write.table(sex, aha3, sep = "\t", quote = F)