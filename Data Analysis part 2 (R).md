### Load packages 
```{r, message=FALSE}
library(ggplot2)
library(ggrepel)
library(tximport)
library(DESeq2)
library(pheatmap)
library(ggvenn)
```
### Read the transcriptome annotation table (from InterProScan)

```{r, message=FALSE}
annotation_data <- read.table("DEG/transcripts_hybrid_annotation.tsv", sep = "\t", quote = "", row.names = NULL, stringsAsFactors = FALSE, fill = TRUE)
#remove orf from the NODe name:
annotation_data$V1 <- sub("_[^_]*$", "", annotation_data$V1)
#selecting PFAM only with highest score for each NODE:
annotation_data$V9 <- gsub("(?<!\\d)\\-(?!\\d)", "0.0", annotation_data$V9, perl = TRUE)
annotation_data$V9 <- as.numeric(annotation_data$V9)
annotation_data <- annotation_data[order(annotation_data$V5),]
annotation_data_filtered <-setDT(annotation_data)[ , .SD[which.max(V9)], by = V1] 
data <- annotation_data_filtered[,-8] # if u want to delete some columns
colnames(data)[1] = "NODE"
colnames(data)[3] = "PF"
colnames(data)[8] = "IPR"
```
### Read row counts from Salmon
```
samples <- list.files(full.names = F, pattern="20_.*x|30_.*x|c_.*x")
files <- file.path(samples, "quant.sf")
names(files) <- str_replace(samples, "~", "") %>% str_replace(".salmon", "")

txi <- tximport(files, type="salmon", tx2gene=data[,c("NODE", "IPR")] , countsFromAbundance="no")

sampletype <- factor(c(rep("20_days",6),rep("30_days",6), rep("no_treatment", 6)))
meta <- data.frame(sampletype, row.names = colnames(txi$counts))

all(colnames(txi$counts) %in% rownames(meta)) # Check that sample names match in both files
```

#### Distribution of variance for row reads

```{r}
data_v <- txi$counts %>% 
  round() %>% 
  data.frame()

mean_counts <- apply(data_v[,1:18], 1, mean)
variance_counts <- apply(data_v[,1:18], 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")
```

# DEG analysis
```{r}
#Create DESeq2Dataset object, run count normalizasion and DE analysis:
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sampletype)
dds <- DESeq(dds)
res <- results(dds)
```








> **_NOTE:_** 
