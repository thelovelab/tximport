context("counts_from_abundance")
test_that("getting counts from abundance works", {

  library(readr)
  dir <- system.file("extdata", package="tximportData")
  samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
  files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
  names(files) <- paste0("sample",1:6)
  tx2gene <- read.csv(file.path(dir, "tx2gene.gencode.v27.csv"))

  txi <- tximport(files, type="salmon", tx2gene=tx2gene)
  txi.S <- tximport(files, type="salmon", tx2gene=tx2gene, 
                      countsFromAbundance="scaledTPM")
  txi.LS <- tximport(files, type="salmon", tx2gene=tx2gene, 
                          countsFromAbundance="lengthScaledTPM")
  
  expect_true(ncol(txi.S$counts) == length(files))
  expect_true(ncol(txi.LS$counts) == length(files))

  # also txOut=TRUE

  txi.tx <- tximport(files, type="salmon", txOut=TRUE,
                     countsFromAbundance="no")
  txi.tx.S <- tximport(files, type="salmon", txOut=TRUE,
                       countsFromAbundance="scaledTPM")
  txi.tx.LS <- tximport(files, type="salmon", txOut=TRUE,
                        countsFromAbundance="lengthScaledTPM")

  # these should not be exactly the same
  # lengthScaledTPM is very close, but adjusted for bias
  expect_true(!all(txi$counts[,1] == txi.S$counts[,1]))
  expect_true(!all(txi$counts[,1] == txi.LS$counts[,1]))

  # what if someone sumToGene() with CFA="no" after it was non-no
  expect_warning({
    txi.sum.S <- summarizeToGene(txi.tx.S, tx2gene=tx2gene,
                                 countsFromAbundance="no")
  }, "incoming counts")
  expect_true(txi.sum.S$countsFromAbundance == "scaledTPM")

  expect_warning({
    txi.sum.LS <- summarizeToGene(txi.tx.LS, tx2gene=tx2gene,
                                  countsFromAbundance="no")
  }, "incoming counts")
  expect_true(txi.sum.LS$countsFromAbundance == "lengthScaledTPM")


  # dtuScaledTPM
  txi.tx.dtu <- tximport(files, type="salmon", tx2gene=tx2gene,
                         txOut=TRUE, countsFromAbundance="dtuScaledTPM")

  ## cors <- sapply(seq_len(nrow(txi.tx.S$counts)), function(i) {
  ##   x <- txi.tx.LS$counts[i,]
  ##   y <- txi.tx.dtu$counts[i,]
  ##   if (all(x==0) | all(y==0)) NA else cor(x,y)
  ## })
  ## hist(cors[cors > .98], col="grey")
  
  # errors for these:
  expect_error(tximport(files, type="salmon", txOut=TRUE, countsFromAbundance="dtuScaledTPM"))
  expect_error(tximport(files, type="salmon", tx2gene=tx2gene, txOut=FALSE, countsFromAbundance="dtuScaledTPM"))
  
})
