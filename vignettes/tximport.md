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
txi <- tximport(files, type="kallisto", gene2tx=gene2tx)
```

```
## reading in files
## 1 2 3 4 5 6 
## transcripts missing genes: 3
## summarizing abundance
## summarizing counts
## summarizing length
```

```r
head(txi$counts)
```

```
##                [,1]      [,2]       [,3]      [,4]      [,5]     [,6]
## A1BG     107.612000 314.47100 109.020000 115.00000  85.87900 76.00740
## A1BG-AS1  82.991700 137.25100 127.188000 142.94000 132.74500 91.08470
## A1CF       9.003612  12.00968   3.005026  15.01005  24.01227 22.01550
## A2M       24.000000   2.00000  20.000000   6.00000  38.00000  8.00000
## A2M-AS1    1.000000   1.00000   1.000000   1.00000   0.00000  0.00000
## A2ML1      3.013540   1.01697   3.049390   2.05004   2.02494  3.04791
```

We can also generate counts from abundances and the average transcript length,
averaged over samples. Using this approach, the counts are not correlated 
with length, and so the length matrix does not need to be provided as an offset
for downstream analysis packages.


```r
txi.cfa <- tximport(files, type="kallisto", gene2tx=gene2tx, countsFromAbundance=TRUE)
```

```
## reading in files
## 1 2 3 4 5 6 
## transcripts missing genes: 3
## summarizing abundance
## summarizing counts
## summarizing length
```

```r
head(txi.cfa$counts)
```

```
##                 [,1]        [,2]        [,3]        [,4]       [,5]
## A1BG     107.6911343 313.5336992 109.1138731 114.7410500  85.901938
## A1BG-AS1  83.0554867 136.9280419 127.3222728 142.3122493 133.091351
## A1CF       9.0112581  12.0087111   3.0102302  14.8310027  24.279929
## A2M       24.0203432   1.9985216  20.0311894   5.9410636  38.329447
## A2M-AS1    0.9993389   0.9962978   0.9996455   0.9939186   0.000000
## A2ML1      3.2057484   1.0798091   3.2459302   2.1598171   2.168484
##               [,6]
## A1BG     76.112186
## A1BG-AS1 91.108860
## A1CF     21.932863
## A2M       7.979380
## A2M-AS1   0.000000
## A2ML1     2.345353
```

We can also avoid gene-level summarization:


```r
txi.txout <- tximport(files, type="kallisto", txOut=TRUE)
```

```
## reading in files
## 1 2 3 4 5 6
```

```r
head(txi.txout$counts)
```

```
##                     [,1]        [,2]        [,3]        [,4]        [,5]
## NR_001526    0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00
## NR_001526_1  0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00
## NR_001526_2  0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00
## NM_130786    1.07612e+02 3.14471e+02 1.09020e+02 1.15000e+02 8.58790e+01
## NR_015380    8.29917e+01 1.37251e+02 1.27188e+02 1.42940e+02 1.32745e+02
## NM_001198818 1.18209e-04 7.53928e-05 4.19520e-05 3.03545e-04 3.05113e-04
##                     [,6]
## NR_001526    0.00000e+00
## NR_001526_1  0.00000e+00
## NR_001526_2  0.00000e+00
## NM_130786    7.60074e+01
## NR_015380    9.10847e+01
## NM_001198818 1.58659e-04
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
## reading in files
## 1 2 3 4 5 6 
## transcripts missing genes: 3
## summarizing abundance
## summarizing counts
## summarizing length
```

```r
head(txi.salmon$counts)
```

```
##                [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
## A1BG     109.473000 317.44800 110.83600 116.35000  87.30210  76.46970
## A1BG-AS1  81.548100 134.82700 136.27800 154.00100 137.40000 101.87300
## A1CF       9.035861  11.05221   5.02241  14.03400  25.36073  25.07424
## A2M       24.000000   2.00000  21.00000   6.00000  38.00000   8.00000
## A2M-AS1    1.000000   1.00000   1.00000   1.00000   0.00000   0.00000
## A2ML1      3.075060   1.03979   4.12350   1.07323   2.13262   6.24507
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
## reading in files
## 1 2 3 4 5 6
```

```r
head(txi.rsem$counts)
```

```
##           [,1]   [,2]  [,3]   [,4]  [,5]  [,6]
## A1BG     94.64 278.03 94.07  96.00 55.00 64.03
## A1BG-AS1 64.28 114.08 98.88 109.05 95.32 73.11
## A1CF      0.00   2.00  1.00   1.00  0.00  1.00
## A2M      24.00   2.00 18.00   4.00 35.00  8.00
## A2M-AS1   1.00   1.00  1.00   0.00  0.00  0.00
## A2ML1     0.84   2.89  0.00   1.00  2.00  3.11
```

## Import with edgeR, DESeq2, limma+voom

TODO
