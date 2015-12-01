<!--
%\VignetteEngine{knitr}
%\VignetteIndexEntry{tximport}
-->

# tximport: import and summarize transcript-level estimates for gene-level analysis

## kallisto

First, read in some kallisto example files:


```r
library(tximportData)
dir <- system.file("extdata", package="tximportData")
list.files(dir)
```

```
## [1] "cufflinks"   "gene2tx.csv" "kallisto"    "rsem"        "salmon"     
## [6] "samples.txt" "tmp"
```

```r
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
files <- file.path(dir,"kallisto", samples$run, "abundance.tsv")
file.exists(files)
```

```
## [1] TRUE TRUE TRUE TRUE TRUE TRUE
```

Transcripts need to be associated with gene IDs for summarization.
If that information is present in the files, we can skip this step.
But for kallisto and Salmon, the files just provide the transcript ID.
So we first make a `data.frame` with two columns: gene ID, transcript ID.
This can be accomplished from a *TxDb* object and the `select` function. 
(TODO: show example)


```r
gene2tx <- read.csv(file.path(dir, "gene2tx.csv"))
```

The *tximport* package has a single function for importing transcript-level estimates.
The `type` argument is used to specify what software was used for estimation
("kallisto", "salmon", "rsem" implemented so far).
A simple list with matrices, "abundance", "counts", and "length", is returned.
The "length" matrix can be used to generate an offset matrix for downstream
gene-level differential analysis of count matrices.


```r
library(tximport)
```

```
## Error in library(tximport): there is no package called 'tximport'
```

```r
txi <- tximport(files, type="kallisto", gene2tx=gene2tx)
```

```
## Error in eval(expr, envir, enclos): could not find function "tximport"
```

```r
names(txi)
```

```
## Error in eval(expr, envir, enclos): object 'txi' not found
```

```r
head(txi$counts)
```

```
## Error in head(txi$counts): object 'txi' not found
```

We can also generate counts from abundances and the average transcript length,
averaged over samples. Using this approach, the counts are not correlated 
with length, and so the length matrix does not need to be provided as an offset
for downstream analysis packages.


```r
txi.cfa <- tximport(files, type="kallisto", gene2tx=gene2tx, countsFromAbundance=TRUE)
```

```
## Error in eval(expr, envir, enclos): could not find function "tximport"
```

```r
head(txi.cfa$counts)
```

```
## Error in head(txi.cfa$counts): object 'txi.cfa' not found
```

We can also avoid gene-level summarization:


```r
txi.txout <- tximport(files, type="kallisto", txOut=TRUE)
```

```
## Error in eval(expr, envir, enclos): could not find function "tximport"
```

```r
head(txi.txout$counts)
```

```
## Error in head(txi.txout$counts): object 'txi.txout' not found
```

## Salmon


```r
files <- file.path(dir,"salmon", samples$run, "quant.sf")
file.exists(files)
```

```
## [1] TRUE TRUE TRUE TRUE TRUE TRUE
```

```r
txi.salmon <- tximport(files, type="salmon", gene2tx=gene2tx)
```

```
## Error in eval(expr, envir, enclos): could not find function "tximport"
```

```r
head(txi.salmon$counts)
```

```
## Error in head(txi.salmon$counts): object 'txi.salmon' not found
```

## RSEM


```r
files <- file.path(dir,"rsem", samples$run, paste0(samples$run, ".genes.results"))
file.exists(files)
```

```
## [1] TRUE TRUE TRUE TRUE TRUE TRUE
```

```r
txi.rsem <- tximport(files, type="rsem")
```

```
## Error in eval(expr, envir, enclos): could not find function "tximport"
```

```r
head(txi.rsem$counts)
```

```
## Error in head(txi.rsem$counts): object 'txi.rsem' not found
```

## Import with edgeR, DESeq2, limma+voom

TODO
