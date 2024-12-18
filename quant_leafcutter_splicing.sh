

conda activate base
module load plink2
module load bedtools2
module load samtools
module load tabix


# Extract junction files
rm -i mappedReads/brbSeqJunc/brbSeqjuncFiles.txt

for bamfile in `ls mappedReads/*.sorted.bam`; do
    echo Converting $bamfile to $bamfile.junc
    sampID=$(basename "${bamfile}.junc")
    /home/c/cs806/regtools/build/regtools junctions extract -a 8 -m 50 -M 500000 -s XS $bamfile -o mappedReads/brbSeqJunc/${sampID}
    echo mappedReads/brbSeqJunc/${sampID} >> mappedReads/brbSeqJunc/brbSeqjuncFiles.txt
done


conda activate base
python /home/c/cs806/leafcutter/clustering/leafcutter_cluster_regtools.py -j mappedReads/brbSeqJunc/brbSeqjuncFiles.txt -m 50 -o mappedReads/brbSeqJunc/juncClus -l 500000

module purge
python /home/c/cs806/leafcutter/scripts/prepare_phenotype_table.py mappedReads/brbSeqJunc/juncClus_perind.counts.gz -p 10

# DS script not running
#/home/c/cs806/leafcutter/scripts/leafcutter_ds.R --num_threads 4 mappedReads/brbSeqJunc/juncClus_perind_numers.counts.gz mappedReads/brbSeqJunc/brbSeqjuncFiles.txt
