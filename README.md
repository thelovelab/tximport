# tximport

Import and summarize transcript-level estimates for gene-level analysis

Description of methods and analysis described in:

* Charlotte Soneson, Michael I. Love, Mark D. Robinson.
[Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences](http://f1000research.com/articles/4-1521/v1),
*F1000Research*, 4:1521, December 2015. doi: 10.12688/f1000research.7563.1

---

Imports transcript-level abundance, estimated counts and 
transcript lengths, and summarizes into matrices for use with downstream
gene-level analysis packages such as edgeR, DESeq2, limma-voom. 
Average transcript length, weighted by 
sample-specific transcript abundance estimates, is provided as a matrix
which can be used as an offset for different expression of 
gene-level counts.

See examples in [vignette](https://github.com/mikelove/tximport/blob/master/vignettes/tximport.md)

Notes:

* tximport does not import the bootstrap estimates from kallisto,
  Salmon, or Sailfish
* Though we provide here functionality for performing gene-level
  differential expression using summarized transcript-level estimates,
  this is does not mean we suggest that users *only* perform gene-level
  analysis. Gene-level differential expression can be complemented
  with transcript- or exon-level analysis. The argument `txOut=TRUE`
  can be used to generate transcript-level matrices.
