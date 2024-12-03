# Remove rRNA in SortMeRNA v4.3.4 

#SLURM job:

´´´
for SAMPLE in 3_Ro8x
do
  sortmerna \
    --fastx True\
    --paired_in True \
    -other /home/kuprinak/sortmerna/run/out/${SAMPLE}_sortmerna\
    --out2 True \
    --idx-dir /home/kuprinak/sortmerna/run/idx \
    --ref /home/kuprinak/RNAseq/fastq/rRNA_databases_v4.3.4/smr_v4.3_default_db.fasta \
    --reads /home/kuprinak/RNAseq/fastq/${SAMPLE}_1.fastq.gz \
    --reads /home/kuprinak/RNAseq/fastq/${SAMPLE}_2.fastq.gz \
    --threads 20 \
    -v
done

´´´
