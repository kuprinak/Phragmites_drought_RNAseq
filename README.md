# _Phragmites_ RNAseq
## Transcriptomic analysis of _Phragmites australis_ gene expression under drought stress

### Used programs:

* SortMeRNA v4.3.4 - filtering out rRNA
* Trimmomatic v0.33 - trimming of reads for adapters, low-quality bases and low-quality reads
* Kraken2 - filtering out conatminating reads
* FastQC v0.12.0 - read quality report
* rnaSPAdes v3.15.4  - transcriptome assembling
* RNAquast v2.2.3 and BUSCO v5.4.4 - transcriptome quality assessment
* KAAS v2.1 and InterProScan v5.61-93.0 - transcriptome annotation
* Salmon v1.10.3 - read quantification
* DESeq2 v1.44.0 - counts normalisation and differential gene expression analysis
* clusterProfiler v4.12.6 - enrichment analysis for differentially expressed genes

### Workflow:
```mermaid
flowchart TB
    A@{shape: procs, label: "Illumina raw reads"} --> B[SortMeRNA];
AA@{shape: procs, label: "Nanopore raw reads"} --> G[rnaSPAdes];
K2 -->T@{shape: cyl, label: "_De novo_ Transcriptome"};
    B --> C[Trimmomatic];
    C --> K[Kraken2];
    K --> E[FastQC];
    K --> F@{shape: procs, label: "Clean reads"};
    A --> G[rnaSPAdes];
    G --> K2[Kraken2];
    T --> BU[BUSCO];
    T --> Q[rnaQUAST];
    T --> I[InterProScan];
    I --> An@{shape: procs, label: "Annotation"};
    KAAS --> An@{shape: procs, label: "Annotation"};
    An --> R;
    T --> KAAS[KAAS];  
    T --> S[Salmon];
    F --> S[Salmon];
    S --> q@{shape: procs, label: "counts"};
    q --> R[DESeq2];
    R --> DEG@{shape: procs, label: "Differentially expressed genes"};
    DEG --> cl[clusterProfiler];
```
