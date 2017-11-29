context("no_tx2gene")
test_that("no tx2gene provided throws error", {

  dir <- system.file("extdata", package="tximportData")
  samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
  files <- file.path(dir,"salmon", samples$run, "quant.sf")
  names(files) <- paste0("sample",1:6)
  expect_error(tximport(files, type="salmon", txOut=FALSE))

})
