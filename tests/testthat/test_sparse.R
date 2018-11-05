context("sparse")
test_that("importing sparsely works", {

  library(readr)
  dir <- system.file("extdata", package="tximportData")
  samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
  files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
  names(files) <- paste0("sample",1:6)

  tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))

  txi0 <- tximport(files, type="salmon", txOut=TRUE)
  txi <- tximport(files, type="salmon", txOut=TRUE, sparse=TRUE)
  idx <- txi0$counts[,1] >= 1
  expect_equal(txi0$counts[idx,1], txi$counts[idx,1])

  txi.cfa0 <- tximport(files, type="salmon", txOut=TRUE, countsFromAbundance="scaledTPM")
  txi.cfa <- tximport(files, type="salmon", txOut=TRUE, countsFromAbundance="scaledTPM", sparse=TRUE)
  idx <- txi0$counts[,1] >= 1
  # test for equality with some tolerance (not exactly equal bc of thresholding for counts < 1)
  expect_equal(txi.cfa0$counts[idx,1], txi.cfa$counts[idx,1], tolerance=.1)  
  
})
