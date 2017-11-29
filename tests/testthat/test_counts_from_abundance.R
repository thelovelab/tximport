context("counts_from_abundance")
test_that("getting counts from abundance works", {

  library(readr)
  dir <- system.file("extdata", package="tximportData")
  samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
  files <- file.path(dir,"salmon", samples$run, "quant.sf")
  names(files) <- paste0("sample",1:6)
  tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))
  
  txi.cfa <- tximport(files, type="salmon", tx2gene=tx2gene, 
                      countsFromAbundance="scaledTPM")
  txi.cfa.len <- tximport(files, type="salmon", tx2gene=tx2gene, 
                          countsFromAbundance="lengthScaledTPM")
  
  expect_true(ncol(txi.cfa$counts) == length(files))
  expect_true(ncol(txi.cfa.len$counts) == length(files))

  # also txOut=TRUE

  txi <- tximport(files, type="salmon", txOut=TRUE,
                  countsFromAbundance="no")
  txi.cfa <- tximport(files, type="salmon", txOut=TRUE,
                      countsFromAbundance="scaledTPM")
  txi.cfa.len <- tximport(files, type="salmon", txOut=TRUE,
                          countsFromAbundance="lengthScaledTPM")

  # these should not be exactly the same
  # lengthScaledTPM is very close, but adjusted for bias
  expect_true(!all(txi$counts[,1] == txi.cfa$counts[,1]))
  expect_true(!all(txi$counts[,1] == txi.cfa.len$counts[,1]))

})
