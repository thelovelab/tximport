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
names(files) <- paste0("sample",1:6)
```

Transcripts need to be associated with gene IDs for summarization.
If that information is present in the files, we can skip this step.
But for kallisto and Salmon, the files just provide the transcript ID.
So we first make a `data.frame` with two columns: gene ID (column 1)
and transcript ID (column 2).
The column names are not relevant but this column order must be used.
This can be accomplished from a *TxDb* object and the `select` function. 
(TODO: show example)


```r
gene2tx <- read.csv(file.path(dir, "gene2tx.csv"))
head(gene2tx)
```

```
##     GENEID       TXNAME
## 1     A1BG    NM_130786
## 2 A1BG-AS1    NR_015380
## 3     A1CF NM_001198818
## 4     A1CF NM_001198819
## 5     A1CF NM_001198820
## 6     A1CF    NM_014576
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
names(txi)
```

```
## [1] "abundance"           "counts"              "length"             
## [4] "countsFromAbundance"
```

```r
head(txi$counts)
```

```
##             sample1   sample2    sample3   sample4   sample5  sample6
## A1BG     107.612000 314.47100 109.020000 115.00000  85.87900 76.00740
## A1BG-AS1  82.991700 137.25100 127.188000 142.94000 132.74500 91.08470
## A1CF       9.003612  12.00968   3.005026  15.01005  24.01227 22.01550
## A2M       24.000000   2.00000  20.000000   6.00000  38.00000  8.00000
## A2M-AS1    1.000000   1.00000   1.000000   1.00000   0.00000  0.00000
## A2ML1      3.013540   1.01697   3.049390   2.05004   2.02494  3.04791
```

We can also generate counts from abundances, scaled to library size
(scaledTPM) or additionally scaled using the average transcript length,
averaged over samples and to library size (lengthScaledTPM). 
Using either of these approaches, the counts are not correlated 
with length, and so the length matrix does not need to be provided as an offset
for downstream analysis packages.


```r
txi.cfa <- tximport(files, type="kallisto", gene2tx=gene2tx, countsFromAbundance="scaledTPM")
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
txi.cfa.len <- tximport(files, type="kallisto", gene2tx=gene2tx, countsFromAbundance="lengthScaledTPM")
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
##              sample1     sample2     sample3     sample4     sample5
## A1BG     103.4635896 270.0794930 109.3824052 102.8084375  91.9913919
## A1BG-AS1  65.6857711  97.0946198 105.0672256 104.9657484 117.3246990
## A1CF       1.5949241   1.9056836   0.5559224   2.4480905   4.7900345
## A2M        8.5356651   0.6367469   7.4271946   1.9689048  15.1819430
## A2M-AS1    0.7701529   0.6884195   0.8038409   0.7143603   0.0000000
## A2ML1      1.3208395   0.3989031   1.3954674   0.8299271   0.9958949
##            sample6
## A1BG     76.457522
## A1BG-AS1 75.339376
## A1CF      4.058900
## A2M       2.964734
## A2M-AS1   0.000000
## A2ML1     1.010386
```

```r
head(txi.cfa.len$counts)
```

```
##              sample1     sample2     sample3     sample4    sample5
## A1BG     107.6911343 313.5336992 109.1138731 114.7410500  85.901938
## A1BG-AS1  83.0554867 136.9280419 127.3222728 142.3122493 133.091351
## A1CF       9.0112581  12.0087111   3.0102302  14.8310027  24.279929
## A2M       24.0203432   1.9985216  20.0311894   5.9410636  38.329447
## A2M-AS1    0.9993389   0.9962978   0.9996455   0.9939186   0.000000
## A2ML1      3.2057484   1.0798091   3.2459302   2.1598171   2.168484
##            sample6
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
##                  sample1     sample2     sample3     sample4     sample5
## NR_001526    0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00
## NR_001526_1  0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00
## NR_001526_2  0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00
## NM_130786    1.07612e+02 3.14471e+02 1.09020e+02 1.15000e+02 8.58790e+01
## NR_015380    8.29917e+01 1.37251e+02 1.27188e+02 1.42940e+02 1.32745e+02
## NM_001198818 1.18209e-04 7.53928e-05 4.19520e-05 3.03545e-04 3.05113e-04
##                  sample6
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
names(files) <- paste0("sample",1:6)
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
##             sample1   sample2   sample3   sample4   sample5   sample6
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
names(files) <- paste0("sample",1:6)
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
##          sample1 sample2 sample3 sample4 sample5 sample6
## A1BG       94.64  278.03   94.07   96.00   55.00   64.03
## A1BG-AS1   64.28  114.08   98.88  109.05   95.32   73.11
## A1CF        0.00    2.00    1.00    1.00    0.00    1.00
## A2M        24.00    2.00   18.00    4.00   35.00    8.00
## A2M-AS1     1.00    1.00    1.00    0.00    0.00    0.00
## A2ML1       0.84    2.89    0.00    1.00    2.00    3.11
```

## Import with edgeR, DESeq2, limma-voom

An example of creating a `DGEList` for use with edgeR:


```r
library(edgeR)
```


```r
cts <- txi$counts
normMat <- txi$length
normMat <- normMat / exp(rowMeans(log(normMat)))
library(edgeR)
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y$offset <- t(t(log(normMat)) + o)
# y is now ready for estimate dispersion functions
# see edgeR User's Guide
```

An example of creating a `DESeqDataSet` for use with DESeq2
(requires version >= 1.11.6).


```r
library(DESeq2)
```


```r
sampleTable <- data.frame(condition=factor(rep(c("A","B"),each=3)))
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ condition)
```

```
## using counts and average transcript lengths from tximport
```

```r
# dds is now ready for DESeq()
# see DESeq2 vignette
```

An example for use with limma-voom:


```r
library(limma)
design <- model.matrix(~ condition, data=sampleTable)
v <- voom(y, design)
# v is now ready for lmFit()
# see limma User's Guide
```

## Session info


```r
sessionInfo()
```

```
## R Under development (unstable) (2015-07-02 r68623)
## Platform: x86_64-apple-darwin14.3.0 (64-bit)
## Running under: OS X 10.10.5 (Yosemite)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices datasets  utils    
## [8] methods   base     
## 
## other attached packages:
##  [1] tximport_0.0.7             DESeq2_1.11.7             
##  [3] RcppArmadillo_0.5.600.2.0  Rcpp_0.12.1               
##  [5] SummarizedExperiment_0.3.9 Biobase_2.29.1            
##  [7] GenomicRanges_1.21.29      GenomeInfoDb_1.5.16       
##  [9] IRanges_2.3.22             S4Vectors_0.7.18          
## [11] BiocGenerics_0.15.6        edgeR_3.11.3              
## [13] limma_3.25.16              tximportData_0.1          
## [15] knitr_1.11                 devtools_1.9.1            
## 
## loaded via a namespace (and not attached):
##  [1] genefilter_1.51.1     locfit_1.5-9.1        reshape2_1.4.1       
##  [4] splines_3.3.0         lattice_0.20-33       colorspace_1.2-6     
##  [7] survival_2.38-3       XML_3.98-1.3          foreign_0.8-66       
## [10] DBI_0.3.1             BiocParallel_1.3.52   RColorBrewer_1.1-2   
## [13] lambda.r_1.1.7        plyr_1.8.3            stringr_1.0.0        
## [16] zlibbioc_1.15.0       munsell_0.4.2         gtable_0.1.2         
## [19] futile.logger_1.4.1   codetools_0.2-14      memoise_0.2.1        
## [22] evaluate_0.8          latticeExtra_0.6-26   geneplotter_1.47.0   
## [25] AnnotationDbi_1.31.18 proto_0.3-10          acepack_1.3-3.3      
## [28] xtable_1.7-4          scales_0.3.0          formatR_1.2.1        
## [31] Hmisc_3.17-0          annotate_1.47.4       XVector_0.9.4        
## [34] gridExtra_2.0.0       ggplot2_1.0.1         digest_0.6.8         
## [37] stringi_0.5-5         grid_3.3.0            tools_3.3.0          
## [40] magrittr_1.5          RSQLite_1.0.0         Formula_1.2-1        
## [43] cluster_2.0.3         futile.options_1.0.0  MASS_7.3-44          
## [46] rpart_4.1-10          nnet_7.3-11           compiler_3.3.0       
## [49] git2r_0.11.0
```
