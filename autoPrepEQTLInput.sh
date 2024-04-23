#!/bin/bash

# This script takes two positional arguments
# - path to plink file to be processed
# - name prefix for processed genotype file
# The scripts should be run from the output directory


echo ""
echo ""
echo "Path to the genotype data in plink format to be processed : $1"
echo "Processed output files will have the prefix : $2"
echo ""
echo ""

# SNP dosage
plink --bfile $1 --recode A-transpose --geno 0.1 --mind 0.5 --maf 0.01 --hwe 1e-6 --out temp --allow-extra-chr
# SNP location
cut -f 2,7- temp.traw > genoFiles/$2.Geno.txt
# Which allele is counted in the dataset
awk '{print$2,$1,$4}' temp.traw > genoFiles/$2.Geno_snpsloc.txt
# Clean-up
cut -f 1-6 temp.traw > genoFiles/$2.Geno_counted_allele.txt
rm temp*

# # Prepare files for eigenMT
# for N in {1..22}
# do
# plink --bfile $1 --chr $N --make-bed --out genoFiles/TEMP_chr$N --allow-extra-chr
# # SNP dosage
# plink --bfile genoFiles/TEMP_chr$N --recode A-transpose --out genoFiles/TEMP2_chr$N
# cut -f 2,7- genoFiles/TEMP2_chr$N.traw > eigenMT/inputGeno/${2}_SNP_chr$N.txt
# # SNP location
# awk '{print$2,$1,$4}' genoFiles/TEMP2_chr$N.traw > eigenMT/inputGeno/${2}_snpsloc_chr$N.txt
# 
# done
# 
# rm genoFiles/TEMP*

# Perform PCA and get Freq
plink --bfile $1 --geno 0.1 --mind 0.5 --maf 0.01 --hwe 1e-6 --make-bed --out temp --allow-extra-chr
plink --bfile temp --pca  --out genoFiles/$2 --allow-extra-chr
plink --bfile temp --freq --out genoFiles/$2 --allow-extra-chr
rm temp*

module load R
Rscript brbSeq_scripts/autoPrepGenoPCA.R $2
Rscript brbSeq_scripts/autoPrepSex.R $1 $2


# Normalize Expr and harmonize with samples in Geno
Rscript brbSeq_scripts/autoPrepGenoExprCov.R $2

# Get peer factors of expression data
conda deactivate
conda activate peer
Rscript brbSeq_scripts/autoPeer.R $2
conda deactivate


#plink --bfile /scratch/cellfunc/shared/HUVEC_BackUp_Oct22/HUVEC_genotype/Imputed_Genotypes_Data_HUVEC_Samples_with_RNA-seq_Data --keep /scratch/cellfunc/shared/HUVEC_BackUp_Oct22/HUVEC_genotype/huvec50.txt --allow-extra-chr --make-bed --out /scratch/vasccell/cs806/brbSeq/genoFiles/huvec50Geno
