context("kallisto")
test_that("import kallisto works", {

  library(readr)
  dir <- system.file("extdata", package="tximportData")
  samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
  files <- file.path(dir,"kallisto", samples$run, "abundance.tsv")
  names(files) <- paste0("sample",1:6)
  tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))

  txi <- tximport(files, type="kallisto", tx2gene=tx2gene)
  expect_true(ncol(txi$counts) == length(files))
  
})
