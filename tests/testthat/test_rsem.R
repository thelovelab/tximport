context("rsem")
test_that("import rsem works", {

  library(readr)
  dir <- system.file("extdata", package="tximportData")
  samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
  files <- file.path(dir,"rsem", samples$run, paste0(samples$run, ".genes.results.gz"))
  names(files) <- paste0("sample",1:6)
  expect_message(txi.rsem <- tximport(files, type="rsem", txOut=FALSE), "looks like you")

  files <- file.path(dir,"rsem", samples$run, paste0(samples$run, ".isoforms.results.gz"))
  names(files) <- paste0("sample",1:6)
  txi.rsem <- tximport(files, type="rsem", txIn=TRUE, txOut=TRUE)
  
})
