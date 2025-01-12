openjdk version "23.0.1" 2024-10-15

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
FASTQ_DIR="/home/kuprinak/RNAseq/fastq-trimmomatic/PHRED30_LEN100/Trimmed_SortMeRNA"
KRAKEN_DIR="."
OUTPUT_DIR="krakentools_filter"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Get the list of Kraken2 output files and pick the one corresponding to the current task ID
KRAKEN_FILES=($KRAKEN_DIR/*_kraken2_output.txt)
KRAKEN_OUTPUT=${KRAKEN_FILES[$SLURM_ARRAY_TASK_ID]}

# Extract the base name of the sample
SAMPLE=$(basename "$KRAKEN_OUTPUT" _kraken2_output.txt)
echo "Processing sample: $SAMPLE"
FORWARD_FASTQ="$FASTQ_DIR/${SAMPLE}_forward_paired.fq.gz"
REVERSE_FASTQ="$FASTQ_DIR/${SAMPLE}_reverse_paired.fq.gz"

# Step 1: Extract reads of taxa (from "taxid" list) for filtering:
if [[ -f "$FORWARD_FASTQ" && -f "$REVERSE_FASTQ" ]]; then
    echo "Filtering paired-end reads for $SAMPLE..."
    python3 /home/saenkos/kraken2/KrakenTools/extract_kraken_reads.py \
        -k $KRAKEN_OUTPUT \
        -s1 $FORWARD_FASTQ \
        -s2 $REVERSE_FASTQ \
        --taxid 33154 2 1783272 10239 2157  \
        --exclude \
        --include-children \
        --fastq-output \
        --report /home/saenkos/kraken2/kraken2_results_PLUSPFP/${SAMPLE}_kraken2_report.txt \
        -o $OUTPUT_DIR/${SAMPLE}_forward_clean.fq \
        -o2 $OUTPUT_DIR/${SAMPLE}_reverse_clean.fq

    echo "Filtered reads written to: $OUTPUT_DIR/${SAMPLE}_forward_clean.fq.gz and $OUTPUT_DIR/${SAMPLE}_reverse_clean.fq.gz"
else
    echo "FASTQ files for $SAMPLE not found. Skipping..."
fi

```

### Taxid which were filtered out:

* 33154 - Opisthokonta
* 2 - Bacteria
* 2157 - Archaea
* 10239 - Viruses
* 1783272 - Terrabacteria group 	



## 4. Read quality assessment in FastQC v0.12.0  and MultiQC v1.14 

```{bash}
fastqc -t 96 /home/kuprinak/RNAseq/fastq-kraken/*_clean.fq -o /home/kuprinak/RNAseq/FastQC/sortmerna_trimmomatic_clean/

multiqc /home/kuprinak/RNAseq/FastQC/sortmerna_trimmomatic_clean --interactive 
```

## 5. Transcriptome _De novo_ assembly using Illumina short reads (eight control samples) and Nanopore long reads (one sample) in rnaSPAdes-3.15.4
```{bash}
spades.py --pe1-1 /home/kuprinak/RNAseq/fastq-kraken/c_Hu4x_forward_clean.fq \
          --pe1-2 /home/kuprinak/RNAseq/fastq-kraken/c_Hu4x_reverse_clean.fq\
	  --pe2-1 /home/kuprinak/RNAseq/fastq-kraken/c_Hu8x_forward_clean.fq \
          --pe2-2 /home/kuprinak/RNAseq/fastq-kraken/c_Hu8x_reverse_clean.fq \
	  --pe3-1 /home/kuprinak/RNAseq/fastq-kraken/c_Ro4x_forward_clean.fq \
          --pe3-2 /home/kuprinak/RNAseq/fastq-kraken/c_Ro4x_reverse_clean.fq \
	  --pe4-1 /home/kuprinak/RNAseq/fastq-kraken/c_Ro8x_forward_clean.fq \
          --pe4-2 /home/kuprinak/RNAseq/fastq-kraken/c_Ro8x_reverse_clean.fq \
	  --pe5-1 /home/kuprinak/RNAseq/fastq-kraken/c_Ru4x_forward_clean.fq \
          --pe5-2 /home/kuprinak/RNAseq/fastq-kraken/c_Ru4x_reverse_clean.fq \
	  --pe6-1 /home/kuprinak/RNAseq/fastq-kraken/c_Ru8x_forward_clean.fq \
          --pe6-2 /home/kuprinak/RNAseq/fastq-kraken/c_Ru8x_reverse_clean.fq \
          --nanopore /home/kuprinak/RNAseq/transcriptome/SPAdes/nanopore_hybrid2_all_c_U2T_filtered/Ro4x_nanopore.U2T_clean.fq \
          -o /home/kuprinak/RNAseq/transcriptome/SPAdes/nanopore_hybrid_Ro4x_clean/ \
        --rna -t 2 --cov-cutoff off      
#no more than 2 threads!      
```



pyfasta info for transcriptome:
261.204M basepairs in 199141 sequences

## 6. Assembly quality in rnaQUAST v2.2.3 
```{bash}
python rnaQUAST.py -t 96 \
-o /home/kuprinak/RNAseq/transcriptome/SPAdes/rnaQUAST/ \
--transcripts /home/kuprinak/RNAseq/transcriptome/SPAdes/nanopore_hybrid_Ro4x_clean/Transcriptome_hybrid_clean.fasta
```
Output: 
METRICS/TRANSCRIPTS                                    Transcriptome_hybrid_clean  

 == BASIC TRANSCRIPTS METRICS (calculated without reference genome and gene database) == 
Transcripts                                            199141                      
Transcripts > 500 bp                                   121232                      
Transcripts > 1000 bp                                  92608                       
Average length of assembled transcripts                1311.655                    
Longest transcript                                     14692                       
Total length                                           261204350                   
Transcript N50                                         2263  

## 7. Assembly quality in BUSCO v5.4.4 
```{bash}
busco -i /home/kuprinak/RNAseq/transcriptome/SPAdes/nanopore_hybrid_Ro4x_clean/Transcriptome_hybrid_clean.fasta \
		-o BUSCO \
		-m tran \
		-l poales_odb10 \
		-c 96 \
                --augustus    \
		--long \
		--update-data
```
#### BUSCO Output:

| n     |                                     |  %    |
| :---: | ----------------------------------- | :---: |
| 4443	|Complete BUSCOs                      | 90.7 |
| 808	|Complete and single-copy BUSCOs      | 16.5 |
| 3635	|Complete and duplicated BUSCOs       | 74.2 |
| 82	|Fragmented BUSCOs                    | 1.7  |
| 371	|Missing BUSCOs                       | 7.6  |
| 4896	|Total BUSCO groups searched          | 100.0 |

Dependency: hmmsearch v3.1



## 8. Annotation of transcriptome in InterProScan
```{bash}
# Fasta file of the transcriptome was chunked into 10 pieces:
pyfasta split -n 7 /home/kuprinak/RNAseq/transcriptome/SPAdes/nanopore_hybrid_Ro4x_clean/chunked/Transcriptome_hybrid_clean.fasta

# Each chunk was run in interproscan, maximum cpu - 40, otherwise java jump error, took one day for about 10000 sequences:
/home/kuprinak/my_interproscan/interproscan-5.61-93.0/interproscan.sh -cpu 40 -t n -i /home/kuprinak/RNAseq/transcriptome/SPAdes/nanopore_hybrid_Ro4x_clean/chunked/Transcriptome_hybrid_clean.0.fasta  \
                 -b /home/kuprinak/RNAseq/transcriptome/SPAdes/Annotation/Clean_annotation/

#  Output files were merged:
cat *.tsv > transcripts_annotation.tsv
```

## 9. Read quantification in Salmon 
```{bash}
# Building an index for the transcriptome:
salmon index -t /home/kuprinak/RNAseq/transcriptome/SPAdes/transcriptome/transcripts.fasta -i transcripts_clean_index

# Quantification:
for SAMPLE in  2_Hu4x 2_Hu8x 2_Ru4x 2_Ru8x 2_Ro4x 2_Ro8x \
               c_Hu4x c_Hu8x c_Ru4x c_Ru8x c_Ro4x c_Ro8x \
               3_Hu4x 3_Hu8x 3_Ru4x 3_Ru8x 3_Ro4x 3_Ro8x 

do

salmon quant -i /home/kuprinak/jobs/transcripts_clean_index/ -l ISF \
         -1 /home/kuprinak/RNAseq/fastq-kraken/${SAMPLE}_forward_clean.fq \
         -2 /home/kuprinak/RNAseq/fastq-kraken/${SAMPLE}_reverse_clean.fq \
         -p 96 --gcBias --seqBias --posBias --validateMappings -o /home/kuprinak/RNAseq/salmon/Hybrid/clean/${SAMPLE}_quant
done 

```

> **_NOTE:_** Resulted read counts were analysed in R. See "Data Analysis part 2 (R)"
