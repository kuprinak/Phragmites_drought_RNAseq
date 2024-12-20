# _Phragmites_ RNAseq
## Transcriptomic analysis of _Phragmites australis_ gene expression under drought stress

### Used programs:

* [SortMeRNA](https://academic.oup.com/bioinformatics/article/28/24/3211/246053?login=true) v4.3.4 - filtering out rRNA
* [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) v0.33 - trimming of adapters, low-quality bases and removing of low-quality reads
* [Kraken2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) and [KrakenTools](https://github.com/jenniferlu717/KrakenTools) v1.2 - filtering out contaminating reads
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.12.0 and [MultiQC](https://seqera.io/multiqc/) v1.14 - read quality report
* [rnaSPAdes](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03614-2) v3.15.4  - transcriptome assembling
* [RNAquast](https://academic.oup.com/bioinformatics/article/32/14/2210/1743439) v2.2.3 and [BUSCO](https://busco.ezlab.org/) v5.4.4 - transcriptome quality assessment
* [KAAS](https://www.genome.jp/kegg/kaas/) v2.1 and [InterProScan](https://interproscan-docs.readthedocs.io/en/latest/Introduction.html) v5.61-93.0 - transcriptome annotation
* [Salmon](https://combine-lab.github.io/salmon/) v1.10.3 - read quantification
* [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) v1.44.0 R package - read counts normalisation and differential gene expression analysis
* [clusterProfiler](https://pmc.ncbi.nlm.nih.gov/articles/PMC3339379/) v4.12.6 R package - enrichment analysis for differentially expressed genes

### Workflow:
```mermaid
flowchart TB
    A@{shape: procs, label: "Illumina raw reads"} --> B([SortMeRNA]);
AA@{shape: procs, label: "Nanopore raw reads"} --> G([rnaSPAdes]);
K2 -->T@{shape: cyl, label: "_De novo_ transcriptome assembly"};
    B --> C([Trimmomatic]);
    C --> K([Kraken2 + KrakenTools]);
    K --> F@{shape: procs, label: "Clean reads"};
    A --> G([rnaSPAdes]);
    G --> K2([Kraken2 + KrakenTools]);
    T --> BU([BUSCO]);
    T --> Q([rnaQUAST]);
    T --> I([InterProScan]);
    I --> An@{shape: procs, label: "Annotation"};
    KAAS --> An@{shape: procs, label: "Annotation"};
    An --> R;
    T --> KAAS([KAAS]);  
    T --> S([Salmon]);
    F --> S([Salmon]);
    S --> q@{shape: procs, label: "Read counts"};
    q --> R([DESeq2]);
    R --> DEG@{shape: procs, label: "Differentially expressed genes"};
    DEG --> cl([clusterProfiler]);
    K --> E([FastQC]);
    E --> m([MultiQC]);
```
