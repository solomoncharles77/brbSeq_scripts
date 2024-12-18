
module load plink2
module load bcftools
module load samtools
module load tabix


# QC genotype data
plink2 --bfile genoFiles/BRB-seq_Pilot_Study_All_SNPs_No_Filtering --geno 0.05 --hwe 1e-06 --maf 0.25 --recode vcf --out genoFiles/brbSeqGeno 

bcftools view genoFiles/brbSeqGeno.vcf --threads 6 -Oz -o genoFiles/brbSeqGeno.vcf.gz
bcftools annotate -x INFO genoFiles/brbSeqGeno.vcf.gz | bcftools +fill-tags  > genoFiles/brbSeqGeno_QCed_tem.vcf

bcftools query -l genoFiles/brbSeqGeno.vcf.gz > genoFiles/brbSeqGeno_samples.txt

# rename to match expr
awk -F'_' '{print $0, $2}' genoFiles/brbSeqGeno_samples.txt > genoFiles/brbSeqGeno_samples2.txt
bcftools reheader -s genoFiles/brbSeqGeno_samples2.txt -o genoFiles/brbSeqGeno_QCed_tem2.vcf genoFiles/brbSeqGeno_QCed_tem.vcf
bcftools query -l genoFiles/brbSeqGeno_QCed_tem2.vcf > genoFiles/brbSeqGeno_samples3.txt

bcftools view genoFiles/brbSeqGeno_QCed_tem2.vcf -S genoFiles/brbSeqGeno_samples3.txt | sed 's/chr//g' | bgzip > genoFiles/brbSeqGeno_QCed.vcf.gz
bcftools index -t genoFiles/brbSeqGeno_QCed.vcf.gz
rm genoFiles/brbSeqGeno_QCed_tem.vcf genoFiles/brbSeqGeno_QCed_tem2.vcf

# Get MAF ----------------------
bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%MAF\n' genoFiles/brbSeqGeno_QCed.vcf.gz > genoFiles/brbSeqGeno_QCed_MAF.txt

# Get SNP dosage ----------------------------
plink2 --vcf genoFiles/brbSeqGeno_QCed.vcf.gz --recode A-transpose --out temp --allow-extra-chr
# SNP location
cut -f 2,7- temp.traw > genoFiles/brbSeqGeno_QCed_Geno.txt
# Which allele is counted in the dataset
awk '{print$2,$1,$4}' temp.traw > genoFiles/brbSeqGeno_QCed_Geno_snpsloc.txt
# Clean-up
cut -f 1-6 temp.traw > genoFiles/brbSeqGeno_QCed_Geno_counted_allele.txt
rm temp*
head -n 1 genoFiles/brbSeqGeno_QCed_Geno.txt > genoFiles/brbSeqGeno_QCed_MatrixHeader.txt

# Prep Geno covariates
plink2 --vcf genoFiles/brbSeqGeno_QCed.vcf.gz --pca --out covFiles/brbSeqGeno  
Rscript brbSeq_scripts/autoPrepGenoPCA_brbSeq.R

# Prep gene expression bed file ---------------------------
Rscript brbSeq_scripts/prepGeneExprPheno.R

##########################################################################

# Get pheno PCA for TNFA expression
QTLtools pca --bed phenoFiles/geneExpr_TMM_IVTnorm_TNFA_qtlTools_Ready.bed.gz --scale --center --out covFiles/geneExpr_TMM_IVTnorm_TNFA_qtlTools_Ready
# Combine geno and TNFA expression pheno PCA
Rscript brbSeq_scripts/prep_geneExpr_TMM_IVTnorm_TNFA_Cov.R



for chr in {1..22}
do
    echo "Processing chromosome ${chr}"

    # Process VCF files
    bcftools view genoFiles/brbSeqGeno_QCed.vcf.gz --regions ${chr} -Oz -o genoFiles/chr${chr}.vcf.gz
    bcftools index -t genoFiles/chr${chr}.vcf.gz

    # Process expression files
    zgrep -w "^#chr" phenoFiles/geneExpr_TMM_IVTnorm_TNFA_qtlTools_Ready.bed.gz > phenoFiles/chr${chr}.bed
    zgrep -w ^${chr} phenoFiles/geneExpr_TMM_IVTnorm_TNFA_qtlTools_Ready.bed.gz | sort -nk2 >> phenoFiles/chr${chr}.bed
    cat phenoFiles/chr${chr}.bed | bgzip > phenoFiles/chr${chr}.bed.gz
    tabix -f phenoFiles/chr${chr}.bed.gz
    rm phenoFiles/chr${chr}.bed

    # Run Nominal
    QTLtools cis \
        --vcf genoFiles/chr${chr}.vcf.gz \
        --bed phenoFiles/chr${chr}.bed.gz \
        --cov covFiles/geneExpr_TMM_IVTnorm_TNFA_geno3pc_pheno40pc.txt \
        --window 1000000 \
        --std-err \
        --grp-best \
        --out geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL/chr${chr}_geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL_nominal.txt.gz \
        --nominal 1
      
    # Run Permute
    QTLtools cis \
        --vcf genoFiles/chr${chr}.vcf.gz \
        --bed phenoFiles/chr${chr}.bed.gz \
        --cov covFiles/geneExpr_TMM_IVTnorm_TNFA_geno3pc_pheno40pc.txt \
        --window 1000000 \
        --std-err \
        --grp-best \
        --out geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL/chr${chr}_geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL_permute.txt.gz \
        --permute 1000

done


zcat geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL/*_permute.txt.gz | gzip -c > geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL/geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL_permute_ALL.txt.gz
zcat geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL/*_nominal.txt.gz | gzip -c > geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL/geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL_nominal_ALL.txt.gz
zcat geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL/geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL_nominal_ALL.txt.gz | awk -F'\t' '$14 < 0.05' | gzip -c > geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL/geneExpr_TMM_IVTnorm_TNFA_brbSeq_cisEQTL_nominal_ALL_pval0.05.txt.gz


##############################################################################
###############################################################################

# Get pheno PCA for UNTR expression
QTLtools pca --bed phenoFiles/geneExpr_TMM_IVTnorm_UNTR_qtlTools_Ready.bed.gz --scale --center --out covFiles/geneExpr_TMM_IVTnorm_UNTR_qtlTools_Ready
# Combine geno and UNTR expression pheno PCA
Rscript brbSeq_scripts/prep_geneExpr_TMM_IVTnorm_UNTR_Cov.R



for chr in {1..22}
do
    echo "Processing chromosome ${chr}"

    # Process VCF files
    bcftools view genoFiles/brbSeqGeno_QCed.vcf.gz --regions ${chr} -Oz -o genoFiles/chr${chr}.vcf.gz
    bcftools index -t genoFiles/chr${chr}.vcf.gz

    # Process expression files
    zgrep -w "^#chr" phenoFiles/geneExpr_TMM_IVTnorm_UNTR_qtlTools_Ready.bed.gz > phenoFiles/chr${chr}.bed
    zgrep -w ^${chr} phenoFiles/geneExpr_TMM_IVTnorm_UNTR_qtlTools_Ready.bed.gz | sort -nk2 >> phenoFiles/chr${chr}.bed
    cat phenoFiles/chr${chr}.bed | bgzip > phenoFiles/chr${chr}.bed.gz
    tabix -f phenoFiles/chr${chr}.bed.gz
    rm phenoFiles/chr${chr}.bed

    # Run Nominal
    QTLtools cis \
        --vcf genoFiles/chr${chr}.vcf.gz \
        --bed phenoFiles/chr${chr}.bed.gz \
        --cov covFiles/geneExpr_TMM_IVTnorm_UNTR_geno3pc_pheno40pc.txt \
        --window 1000000 \
        --std-err \
        --grp-best \
        --out geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL/chr${chr}_geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL_nominal.txt.gz \
        --nominal 1
        
    # Run Permute
    QTLtools cis \
        --vcf genoFiles/chr${chr}.vcf.gz \
        --bed phenoFiles/chr${chr}.bed.gz \
        --cov covFiles/geneExpr_TMM_IVTnorm_UNTR_geno3pc_pheno40pc.txt \
        --window 1000000 \
        --std-err \
        --grp-best \
        --out geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL/chr${chr}_geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL_permute.txt.gz \
        --permute 1000

done

zcat geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL/*_permute.txt.gz | gzip -c > geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL/geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL_permute_ALL.txt.gz
zcat geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL/*_nominal.txt.gz | gzip -c > geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL/geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL_nominal_ALL.txt.gz
zcat geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL/geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL_nominal_ALL.txt.gz | awk -F'\t' '$14 < 0.05' | gzip -c > geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL/geneExpr_TMM_IVTnorm_UNTR_brbSeq_cisEQTL_nominal_ALL_pval0.05.txt.gz
