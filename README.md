# _Phragmites_ RNAseq
## Transcriptomic analysis of _Phragmites australis_ gene expression under drought stress

### Workflow:

```mermaid
flowchart TD
    A{Raw reads} -->B[SortMeRNA];
    B --> C[Trimmomatic];
    C --> D[Kraken2];
    D --> E[FastQC];
    D --> F{Clean reads};
    F --> G[rnaSPAdes];
```
