# _Phragmites_ RNAseq
## Transcriptomic analysis of _Phragmites australis_ gene expression under drought stress

### Workflow:

```mermaid
flowchart TD
    A[Raw reads]-->B[SortMeRNA];
    B[SortMeRNA]-->C[Trimmomatic];
    D[Trimmomatic]-->D[Kraken2];
    D[Kraken2]->E[FastQC];
```
