<!--
%\VignetteEngine{knitr}
%\VignetteIndexEntry{tximport}
-->

# tximport vignette

Import and summarize transcript-level abundance estimates for
gene-level analysis.



## kallisto

We start with some kallisto TSV files containing transcript
abundance estimates. The pipeline will be nearly
identical for other quantitation tools, which are shown below.
First, read in some kallisto example files:


```r
library(tximportData)
dir <- system.file("extdata", package = "tximportData")
list.files(dir)
```

```
## [1] "cufflinks"            "kallisto"             "rsem"                
## [4] "sailfish"             "salmon"               "samples_extended.txt"
## [7] "samples.txt"          "tmp"                  "tx2gene.csv"
```


```r
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
files <- file.path(dir, "kallisto", samples$run, "abundance.tsv")
names(files) <- paste0("sample", 1:6)
```

Transcripts need to be associated with gene IDs for summarization.
If that information is present in the files, we can skip this step.
But for kallisto, Salmon and Sailfish, the files only provide the transcript ID.
We first make a `data.frame` with two columns: transcript ID (column
1) and gene ID (column 1).
The column names do not matter but this column order must be used.
The transcript ID must be the same one used in the abundance/quantification
files. There is no restriction on the gene ID.

Creating this `data.frame` can be accomplished from a *TxDb* object
and the `select` function from the *AnnotationDbi* package. The
following code is not evaluated, but could be used to construct
such a table:


```r
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]  # tx ID, then gene ID
```

Here we read in a pre-constructed `tx2gene` table:


```r
tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))
head(tx2gene)
```

```
##         TXNAME   GENEID
## 1    NM_130786     A1BG
## 2    NR_015380 A1BG-AS1
## 3 NM_001198818     A1CF
## 4 NM_001198819     A1CF
## 5 NM_001198820     A1CF
## 6    NM_014576     A1CF
```

The *tximport* package has a single function for importing transcript-level estimates.
The `type` argument is used to specify what software was used for estimation
("kallisto", "salmon", "rsem" implemented so far).
A simple list with matrices, "abundance", "counts", and "length", is returned.
The "length" matrix can be used to generate an offset matrix for downstream
gene-level differential analysis of count matrices.

While *tximport* works without any dependencies, it is much faster to
read in files using the *readr* package (version >= 0.2.2).
To use this, we pass the `read_tsv` function to `tximport`.


```r
library(tximport)
library(readr)
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, reader = read_tsv)
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
## A1BG     108.581000 314.42400 110.450000 116.00000  85.80300 75.91360
## A1BG-AS1  86.163600 140.10700 129.994000 146.40800 136.92800 97.29540
## A1CF       9.003863  12.01096   3.005232  15.01082  24.01285 22.01611
## A2M       24.000000   2.00000  21.000000   6.00000  38.00000  8.00000
## A2M-AS1    1.000000   1.00000   1.000000   1.00000   0.00000  0.00000
## A2ML1      3.012760   1.01650   3.049480   2.04965   2.02477  3.04483
```

We could alternatively generate counts from abundances, 
using the argument `countsFromAbundance`,
scaled to library size (`"scaledTPM"`) or additionally scaled 
using the average transcript length,
averaged over samples and to library size (`"lengthScaledTPM"`). 
Using either of these approaches, the counts are not correlated 
with length, and so the length matrix does not need to be provided as an offset
for downstream analysis packages.

We can avoid gene-level summarization by setting `txOut=TRUE`.


```r
txi.tx <- tximport(files, type = "kallisto", txOut = TRUE, tx2gene = tx2gene, 
    reader = read_tsv)
```

These can then be summarized after the fact using the function
`summarizeToGene`, which gives the same matrices as using
`txOut=FALSE` in the first `tximport` call.


```r
txi.sum <- summarizeToGene(txi.tx, tx2gene)
all.equal(txi$counts, txi.sum$counts)
```

```
## [1] TRUE
```

## Salmon / Sailfish

Salmon or Sailfish `quant.sf` files can be imported by setting type to
`"salmon"` or `"sailfish"`.


```r
files <- file.path(dir, "salmon", samples$run, "quant.sf")
names(files) <- paste0("sample", 1:6)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, reader = read_tsv)
head(txi.salmon$counts)
```

```
##             sample1   sample2    sample3   sample4   sample5   sample6
## A1BG     109.232000 316.22400 110.638000 116.00000  86.38430  76.91630
## A1BG-AS1  83.969700 138.44900 119.274000 151.08300 123.98500 103.25100
## A1CF       9.030691  10.01847   5.019242  13.01820  25.21914  25.07356
## A2M       24.000000   2.00000  21.000000   6.00000  38.00000   8.00000
## A2M-AS1    1.000000   1.00000   1.000000   1.00000   0.00000   0.00000
## A2ML1      3.047950   1.02987   4.076160   1.04945   3.07761   5.12409
```

## RSEM


```r
files <- file.path(dir, "rsem", samples$run, paste0(samples$run, ".genes.results"))
names(files) <- paste0("sample", 1:6)
txi.rsem <- tximport(files, type = "rsem", reader = read_tsv)
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

**Note**: there are two suggested ways of importing estimates
for use with gene-level differential expression methods. The first method,
which we show below for *edgeR* and for *DESeq2*, is to use the
estimated counts from the quantification tools, and additionally to use the
transcript-level abundance estimates to calculate an offset that
corrects for changes to the average transcript length across samples.
The code examples below accomplish these steps for you.
The second method is to use the *tximport* argument
`countsFromAbundance="lengthScaledTPM"` or `"scaledTPM"`,
and then to use the count matrix `txi$counts` directly
as you would a regular count matrix with these software.

An example of creating a `DGEList` for use with edgeR:


```r
library(edgeR)
```


```r
cts <- txi$counts
normMat <- txi$length
normMat <- normMat/exp(rowMeans(log(normMat)))
library(edgeR)
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y$offset <- t(t(log(normMat)) + o)
# y is now ready for estimate dispersion functions see edgeR User's Guide
```

An example of creating a `DESeqDataSet` for use with DESeq2
(requires R-devel, Bioconductor 3.3 and DESeq2 version >= 1.11.6, or
you can source the function from the 
[development branch](https://github.com/Bioconductor-mirror/DESeq2/blob/master/R/AllClasses.R#L318-L333)).


```r
library(DESeq2)
```

The user should make sure the rownames of `sampleTable` align with the
colnames of `txi$counts`, if there are colnames. The best practice is
to read `sampleTable` from a CSV file, and to construct `files` from a
column of `sampleTable` before calling `tximport`.


```r
sampleTable <- data.frame(condition = factor(rep(c("A", "B"), each = 3)))
rownames(sampleTable) <- colnames(txi$counts)
```


```r
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
# dds is now ready for DESeq() see DESeq2 vignette
```

An example for use with limma-voom. At the moment, limma-voom does not
use the offset matrix stored in `y$offset`, so we recommend using
the scaled counts generated from abundances, either `"scaledTPM"`
or `"lengthScaledTPM"`:


```r
files <- file.path(dir, "kallisto", samples$run, "abundance.tsv")
names(files) <- paste0("sample", 1:6)
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, reader = read_tsv, 
    countsFromAbundance = "lengthScaledTPM")
library(limma)
y <- DGEList(txi$counts)
y <- calcNormFactors(y)
design <- model.matrix(~condition, data = sampleTable)
v <- voom(y, design)
# v is now ready for lmFit() see limma User's Guide
```

## Session info


```r
sessionInfo()
```

```
## R Under development (unstable) (2015-11-08 r69614)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 15.10
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices datasets  utils    
## [8] methods   base     
## 
## other attached packages:
##  [1] DESeq2_1.11.17              SummarizedExperiment_1.1.18
##  [3] Biobase_2.31.3              GenomicRanges_1.23.13      
##  [5] GenomeInfoDb_1.7.6          IRanges_2.5.24             
##  [7] S4Vectors_0.9.26            BiocGenerics_0.17.3        
##  [9] edgeR_3.13.4                limma_3.27.11              
## [11] readr_0.2.2                 tximport_0.0.19            
## [13] tximportData_0.1.2          knitr_1.12.3               
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.3          compiler_3.3.0       formatR_1.2.1       
##  [4] RColorBrewer_1.1-2   plyr_1.8.3           XVector_0.11.4      
##  [7] tools_3.3.0          zlibbioc_1.17.0      rpart_4.1-10        
## [10] RSQLite_1.0.0        annotate_1.49.0      evaluate_0.8        
## [13] gtable_0.1.2         lattice_0.20-33      DBI_0.3.1           
## [16] gridExtra_2.0.0      genefilter_1.53.1    stringr_1.0.0       
## [19] cluster_2.0.3        locfit_1.5-9.1       grid_3.3.0          
## [22] nnet_7.3-11          AnnotationDbi_1.33.7 XML_3.98-1.3        
## [25] survival_2.38-3      BiocParallel_1.5.16  foreign_0.8-66      
## [28] latticeExtra_0.6-26  Formula_1.2-1        geneplotter_1.49.0  
## [31] ggplot2_2.0.0        magrittr_1.5         Hmisc_3.17-1        
## [34] scales_0.3.0         splines_3.3.0        xtable_1.8-0        
## [37] colorspace_1.2-6     stringi_1.0-1        acepack_1.3-3.3     
## [40] munsell_0.4.2        markdown_0.7.7
```
