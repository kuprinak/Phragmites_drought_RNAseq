# _Phragmites_ RNAseq
## Transcriptomic analysis of _Phragmites australis_ gene expression under drought stress

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
    T --> KEGG[KEGG];  
    T --> S[Salmon];
    F --> S[Salmon];
    S --> q@{shape: procs, label: "quants"};
    q --> R[DESeq2];
```
