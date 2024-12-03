# _Phragmites_ RNAseq
## Transcriptomic analysis of _Phragmites australis_ gene expression under drought stress

### Workflow:

```mermaid
graph TD;
    Raw reads-->SortMeRNA;
    SortMeRNA-->Trimmomatic;
    Trimmomatic-->Kraken2;
    Kraken2->FastQC;
```
