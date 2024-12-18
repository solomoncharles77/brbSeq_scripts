

conda activate base
module load plink2
module load bedtools2
module load samtools
module load tabix


# Extract junction files
rm -i testSeq/mappedReads/brbSeqJunc/brbSeqjuncFiles.txt

for bamfile in `ls testSeq/mappedReads/*.sorted.bam`; do
    echo Converting $bamfile to $bamfile.junc
    sampID=$(basename "${bamfile}.junc")
    /home/c/cs806/regtools/build/regtools junctions extract -a 8 -m 50 -M 500000 -s XS $bamfile -o testSeq/mappedReads/brbSeqJunc/${sampID}
    echo testSeq/mappedReads/brbSeqJunc/${sampID} >> testSeq/mappedReads/brbSeqJunc/brbSeqjuncFiles.txt
done


conda activate base
python /home/c/cs806/leafcutter/clustering/leafcutter_cluster_regtools.py -j testSeq/mappedReads/brbSeqJunc/brbSeqjuncFiles.txt -m 50 -o testSeq/mappedReads/brbSeqJunc/juncClus -l 500000

module purge
python /home/c/cs806/leafcutter/scripts/prepare_phenotype_table.py testSeq/mappedReads/brbSeqJunc/juncClus_perind.counts.gz -p 10

############ Split junction files and prep clusters separately for each group
awk '/TNFa/ {print > "testSeq/mappedReads/brbSeqJunc/brbSeqjuncFiles_TNFa.txt"} /untr/ {print > "testSeq/mappedReads/brbSeqJunc/brbSeqjuncFiles_UNTR.txt"}' testSeq/mappedReads/brbSeqJunc/brbSeqjuncFiles.txt

### TNFA
python /home/c/cs806/leafcutter/clustering/leafcutter_cluster_regtools.py -j testSeq/mappedReads/brbSeqJunc/brbSeqjuncFiles_TNFa.txt -m 50 -o testSeq/mappedReads/brbSeqJunc/TNFA_juncClus -l 500000
python /home/c/cs806/leafcutter/scripts/prepare_phenotype_table.py testSeq/mappedReads/brbSeqJunc/TNFA_juncClus_perind.counts.gz -p 10

### UNTR
python /home/c/cs806/leafcutter/clustering/leafcutter_cluster_regtools.py -j testSeq/mappedReads/brbSeqJunc/brbSeqjuncFiles_UNTR.txt -m 50 -o testSeq/mappedReads/brbSeqJunc/UNTR_juncClus -l 500000
python /home/c/cs806/leafcutter/scripts/prepare_phenotype_table.py testSeq/mappedReads/brbSeqJunc/UNTR_juncClus_perind.counts.gz -p 10

mv *.sorted.gz testSeq/mappedReads/brbSeqJunc/

# Prep intron splicing bed file #########################
Rscript brbSeq_scripts/prep_BRBseq_leafcutterQuant.R

# Prep coveriates -------------------
# Get pheno PCA TNFA
QTLtools pca --bed phenoFiles/LeafcutterSplicing_TNFA_qtlTools_Ready.bed.gz --scale --center --out covFiles/LeafcutterSplicing_TNFA_qtlTools_Ready
# Combine geno and pheno PCA
Rscript brbSeq_scripts/prep_BRBseq_LeafcutterSplicing_TNFA_Cov.R

# Get pheno PCA UNTR
QTLtools pca --bed phenoFiles/LeafcutterSplicing_UNTR_qtlTools_Ready.bed.gz --scale --center --out covFiles/LeafcutterSplicing_UNTR_qtlTools_Ready
# Combine geno and pheno PCA
Rscript brbSeq_scripts/prep_BRBseq_LeafcutterSplicing_UNTR_Cov.R

####################
for chr in {1..22}
do
    echo "Processing chromosome ${chr}"

    # Process VCF files
    bcftools view genoFiles/brbSeqGeno_QCed.vcf.gz --regions ${chr} -Oz -o genoFiles/chr${chr}.vcf.gz
    bcftools index -t genoFiles/chr${chr}.vcf.gz

    # Process expression files
    zgrep -w "^#chr" phenoFiles/LeafcutterSplicing_TNFA_qtlTools_Ready.bed.gz > phenoFiles/chr${chr}.bed
    zgrep -w ^${chr} phenoFiles/LeafcutterSplicing_TNFA_qtlTools_Ready.bed.gz | sort -nk2 >> phenoFiles/chr${chr}.bed
    cat phenoFiles/chr${chr}.bed | bgzip > phenoFiles/chr${chr}.bed.gz
    tabix -f phenoFiles/chr${chr}.bed.gz
    rm phenoFiles/chr${chr}.bed

    # Run QTLtools Nominal
    QTLtools cis \
        --vcf genoFiles/chr${chr}.vcf.gz \
        --bed phenoFiles/chr${chr}.bed.gz \
        --cov covFiles/leafcutterSplicing_TNFA_geno3pc_pheno50pc.txt \
        --window 1000000 \
        --std-err \
        --grp-best \
        --out leafcutterSplicing_TNFA_cisEQTL/chr${chr}_leafcutterSplicing_TNFA_cisEQTL_nominal.txt.gz \
        --nominal 1
        
    # Run QTLtools Permute
    QTLtools cis \
        --vcf genoFiles/chr${chr}.vcf.gz \
        --bed phenoFiles/chr${chr}.bed.gz \
        --cov covFiles/leafcutterSplicing_TNFA_geno3pc_pheno50pc.txt \
        --window 1000000 \
        --std-err \
        --grp-best \
        --out leafcutterSplicing_TNFA_cisEQTL/chr${chr}_leafcutterSplicing_TNFA_cisEQTL_permute.txt.gz \
        --permute 1000

done

zcat leafcutterSplicing_TNFA_cisEQTL/*_permute.txt.gz | gzip -c > leafcutterSplicing_TNFA_cisEQTL/leafcutterSplicing_TNFA_cisEQTL_permute_ALL.txt.gz
zcat leafcutterSplicing_TNFA_cisEQTL/*_nominal.txt.gz | gzip -c > leafcutterSplicing_TNFA_cisEQTL/leafcutterSplicing_TNFA_cisEQTL_nominal_ALL.txt.gz
zcat leafcutterSplicing_TNFA_cisEQTL/leafcutterSplicing_TNFA_cisEQTL_nominal_ALL.txt.gz | awk -F'\t' '$14 < 0.05'| gzip -c > leafcutterSplicing_TNFA_cisEQTL/leafcutterSplicing_TNFA_cisEQTL_nominal_ALL_pval0.05.txt.gz

####################
for chr in {1..22}
do
    echo "Processing chromosome ${chr}"

    # Process VCF files
    bcftools view genoFiles/brbSeqGeno_QCed.vcf.gz --regions ${chr} -Oz -o genoFiles/chr${chr}.vcf.gz
    bcftools index -t genoFiles/chr${chr}.vcf.gz

    # Process expression files
    zgrep -w "^#chr" phenoFiles/LeafcutterSplicing_UNTR_qtlTools_Ready.bed.gz > phenoFiles/chr${chr}.bed
    zgrep -w ^${chr} phenoFiles/LeafcutterSplicing_UNTR_qtlTools_Ready.bed.gz | sort -nk2 >> phenoFiles/chr${chr}.bed
    cat phenoFiles/chr${chr}.bed | bgzip > phenoFiles/chr${chr}.bed.gz
    tabix -f phenoFiles/chr${chr}.bed.gz
    rm phenoFiles/chr${chr}.bed

    # Run QTLtools Nominal
    QTLtools cis \
        --vcf genoFiles/chr${chr}.vcf.gz \
        --bed phenoFiles/chr${chr}.bed.gz \
        --cov covFiles/leafcutterSplicing_UNTR_geno3pc_pheno50pc.txt \
        --window 1000000 \
        --std-err \
        --grp-best \
        --out leafcutterSplicing_UNTR_cisEQTL/chr${chr}_leafcutterSplicing_UNTR_cisEQTL_nominal.txt.gz \
        --nominal 1
        
    # Run QTLtools Permute
    QTLtools cis \
        --vcf genoFiles/chr${chr}.vcf.gz \
        --bed phenoFiles/chr${chr}.bed.gz \
        --cov covFiles/leafcutterSplicing_UNTR_geno3pc_pheno50pc.txt \
        --window 1000000 \
        --std-err \
        --grp-best \
        --out leafcutterSplicing_UNTR_cisEQTL/chr${chr}_leafcutterSplicing_UNTR_cisEQTL_permute.txt.gz \
        --permute 1000

done

zcat leafcutterSplicing_UNTR_cisEQTL/*_permute.txt.gz | gzip -c > leafcutterSplicing_UNTR_cisEQTL/leafcutterSplicing_UNTR_cisEQTL_permute_ALL.txt.gz
zcat leafcutterSplicing_UNTR_cisEQTL/*_nominal.txt.gz | gzip -c > leafcutterSplicing_UNTR_cisEQTL/leafcutterSplicing_UNTR_cisEQTL_nominal_ALL.txt.gz
zcat leafcutterSplicing_UNTR_cisEQTL/leafcutterSplicing_UNTR_cisEQTL_nominal_ALL.txt.gz | awk -F'\t' '$14 < 0.05' | gzip -c > leafcutterSplicing_UNTR_cisEQTL/leafcutterSplicing_UNTR_cisEQTL_nominal_ALL_pval0.05.txt.gz
