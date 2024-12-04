## 1. Remove rRNA in SortMeRNA v4.3.4 

### Indexing of the reference database (rRNA_databases_v4.3.4)
```{bash}
sortmerna \
  --reads /home/kuprinak/RNAseq/fastq/2_Hu4x_1.fastq.gz \
  --reads /home/kuprinak/RNAseq/fastq/2_Hu4x_2.fastq.gz \
  --ref /home/kuprinak/RNAseq/fastq/rRNA_databases_v4.3.4/smr_v4.3_default_db.fasta \
  -other /home/kuprinak/sortmerna/run/out/Clean \
  -aligned /home/kuprinak/sortmerna/run/out/With_rRNA \
  -fastx \
  -paired-in \
  -threads 96

#after this run, all files from the folder sortmerna/run/kvbd have to be removed
```
### Remove rRNA for all samples 
```{bash}
for SAMPLE in  2_Hu4x 2_Hu8x 2_Ru4x 2_Ru8x 2_Ro4x 2_Ro8x \
               c_Hu4x c_Hu8x c_Ru4x c_Ru8x c_Ro4x c_Ro8x \
               3_Hu4x 3_Hu8x 3_Ru4x 3_Ru8x 3_Ro4x 3_Ro8x 
do
  sortmerna \
    --fastx True\
    --paired_in True \
    --other /home/kuprinak/sortmerna/run/out/${SAMPLE}_sortmerna\
    --out2 True \
    --idx-dir /home/kuprinak/sortmerna/run/idx \
    --ref /home/kuprinak/RNAseq/fastq/rRNA_databases_v4.3.4/smr_v4.3_default_db.fasta \
    --reads /home/kuprinak/RNAseq/fastq/${SAMPLE}_1.fastq.gz \
    --reads /home/kuprinak/RNAseq/fastq/${SAMPLE}_2.fastq.gz \
    --threads 20 \
    -v
done

```
## 2. Adapters and quality trimming in Trimmomatic v0.33 
```{bash}
for SAMPLE in  2_Hu4x 2_Hu8x 2_Ru4x 2_Ru8x 2_Ro4x 2_Ro8x \
              c_Hu4x c_Hu8x c_Ru4x c_Ru8x c_Ro4x c_Ro8x \
              3_Hu4x 3_Hu8x 3_Ru4x 3_Ru8x 3_Ro4x 3_Ro8x 
do
trimmomatic PE -threads 4 /home/kuprinak/RNAseq/fastq-sortmerna/${SAMPLE}_sortmerna_fwd.fq.gz /home/kuprinak/RNAseq/fastq-sortmerna/${SAMPLE}_sortmerna_rev.fq.gz \
/home/kuprinak/RNAseq/fastq-trimmomatic/${SAMPLE}_forward_paired.fq.gz \
/home/kuprinak/RNAseq/fastq-trimmomatic/${SAMPLE}_forward_unpaired.fq.gz \
/home/kuprinak/RNAseq/fastq-trimmomatic/${SAMPLE}_reverse_paired.fq.gz \
/home/kuprinak/RNAseq/fastq-trimmomatic/${SAMPLE}_reverse_unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:5:30 MINLEN:100  #cut adapters for illumina, QC with Phred 33 lower than 30, min length 100
done

```

## 3. Contamination removal in Kraken2
```{bash}
Wait for Stepan
```

## 4. Read quality assessment in FastQC v0.12.0  and MultiQC v1.14 

```{bash}
fastqc -t 96 /home/kuprinak/RNAseq/fastq-kraken/*_clean.fq -o /home/kuprinak/RNAseq/FastQC/sortmerna_trimmomatic_clean/

multiqc /home/kuprinak/RNAseq/FastQC/sortmerna_trimmomatic_clean --interactive 
```

## 5. Transcriptome _De novo_ assembly using Illumina short reads (eight control samples) and Nanopore long reads (one sample) in rnaSPAdes-3.15.4
```{bash}
spades.py --pe1-1 /home/kuprinak/RNAseq/fastq/c_Hu4x_1.fastq.gz \
          --pe1-2 /home/kuprinak/RNAseq/fastq/c_Hu4x_2.fastq.gz \
	        --pe2-1 /home/kuprinak/RNAseq/fastq/c_Hu8x_1.fastq.gz \
          --pe2-2 /home/kuprinak/RNAseq/fastq/c_Hu8x_2.fastq.gz \
	        --pe3-1 /home/kuprinak/RNAseq/fastq/c_Ro4x_1.fastq.gz \
          --pe3-2 /home/kuprinak/RNAseq/fastq/c_Ro4x_2.fastq.gz \
	        --pe4-1 /home/kuprinak/RNAseq/fastq/c_Ro8x_1.fastq.gz \
          --pe4-2 /home/kuprinak/RNAseq/fastq/c_Ro8x_2.fastq.gz \
	        --pe5-1 /home/kuprinak/RNAseq/fastq/c_Ru4x_1.fastq.gz \
          --pe5-2 /home/kuprinak/RNAseq/fastq/c_Ru4x_2.fastq.gz \
	        --pe6-1 /home/kuprinak/RNAseq/fastq/c_Ru8x_1.fastq.gz \
          --pe6-2 /home/kuprinak/RNAseq/fastq/c_Ru8x_2.fastq.gz \
          --nanopore /home/kuprinak/RNAseq/nanopore/fastq_Ro4x_nanopore.fq.gz\
          -o /home/kuprinak/RNAseq/transcriptome/SPAdes/transcriptome\
          --rna -t 2 --cov-cutoff off
#no more than 2 threads!      
```

## 6. Assembly quality in rnaQUAST v2.2.3 
```{bash}
python rnaQUAST.py -t 96 \
-o /home/kuprinak/RNAseq/transcriptome/SPAdes/rnaQUAST/ \
--transcripts /home/kuprinak/RNAseq/transcriptome/SPAdes/transcriptome/transcripts.fasta 
```
## 7. Assembly quality in BUSCO v5.4.4 
```{bash}
busco -i /home/kuprinak/RNAseq/transcriptome/SPAdes/transcriptome/transcripts.fasta  \
		-o BUSCO \
		-m tran \
		-l poales_odb10 \
		-c 96 \
                --augustus    \
		--long \
		--update-data
```
### BUSCO Output:




## 8. Annotation of transcriptome in InterProScan
```{bash}
# Fasta file of the transcriptome was chunked into 10 pieces:
pyfasta split -n 10 /home/kuprinak/RNAseq/transcriptome/SPAdes/transcriptome/transcripts.fasta

# Each chunk was run in interproscan, maximum cpu - 40, otherwise java jump error, took one day for about 10000 sequences:
/home/kuprinak/my_interproscan/interproscan-5.61-93.0/interproscan.sh -cpu 40 -t n -i /home/kuprinak/RNAseq/transcriptome/SPAdes/transcriptome/transcripts_1.fasta  \
                 -b /home/kuprinak/RNAseq/transcriptome/SPAdes/Annotation/

#  Output files were merged:
cat *.tsv > transcripts_annotation.tsv
```

## 9. Read quantification in Salmon 
```{bash}
# Building an index for the transcriptome:
salmon index -t /home/kuprinak/RNAseq/transcriptome/SPAdes/transcriptome/transcripts.fasta -i transcripts_index

# Quantification:
for SAMPLE in  2_Hu4x 2_Hu8x 2_Ru4x 2_Ru8x 2_Ro4x 2_Ro8x \
               c_Hu4x c_Hu8x c_Ru4x c_Ru8x c_Ro4x c_Ro8x \
               3_Hu4x 3_Hu8x 3_Ru4x 3_Ru8x 3_Ro4x 3_Ro8x 

do

salmon quant -i /home/kuprinak/jobs/transcripts_clean_index/ -l ISF \
         -1 /home/kuprinak/RNAseq/fastq-kraken/${SAMPLE}_forward_clean.fq \
         -2 /home/kuprinak/RNAseq/fastq-kraken/${SAMPLE}_reverse_clean.fq \
         -p 96 --gcBias --seqBias --posBias --validateMappings -o /home/kuprinak/RNAseq/salmon/${SAMPLE}_quant
done 

```

> **_NOTE:_** Resulted read counts were analysed in R. See "Data Analysis part 2 (R)"
