context("h5")
test_that("kallisto HDF5 import works", {

  library(readr)
  dir <- system.file("extdata", package="tximportData")
  samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
  files <- file.path(dir,"kallisto_boot", samples$run, "abundance.h5")
  names(files) <- paste0("sample",1:6)

  txi <- tximport(files, type="kallisto", txOut=TRUE)
  expect_true("infReps" %in% names(txi))
  txi <- tximport(files, type="kallisto", txOut=TRUE, varReduce=TRUE)
  expect_true("variance" %in% names(txi))
  txi <- tximport(files, type="kallisto", txOut=TRUE, dropInfReps=TRUE)
  expect_true(!any(c("infReps","variance") %in% names(txi)))
  
})
