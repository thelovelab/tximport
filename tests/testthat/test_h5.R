context("h5")
test_that("kallisto HDF5 works", {

  library(readr)
  dir <- system.file("extdata", package="tximportData")
  samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
  files <- file.path(dir,"kallisto_h5", samples$run, "abundance.h5")
  names(files) <- paste0("sample",1:6)
  tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))
  txi <- tximport(files, type="kallisto.h5", txOut=TRUE)
  
})
