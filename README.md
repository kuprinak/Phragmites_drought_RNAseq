# _Phragmites_ RNAseq
## Transcriptomic analysis of _Phragmites australis_ gene expression under drought stress

### Workflow:

```mermaid
flowchart TB
    A{shape: procs, label: "Illumina raw reads"} --> B[SortMeRNA];
    B --> C[Trimmomatic];
    C --> D[Kraken2];
    D --> E[FastQC];
    D --> F{Clean reads};
    A --> G[rnaSPAdes];
    AA{Nanopore Raw reads} --> G[rnaSPAdes];
    G --> T{_De novo_ Transcriptome};
    
```
