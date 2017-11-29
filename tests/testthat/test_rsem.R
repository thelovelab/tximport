context("rsem")
test_that("import rsem works", {
  library(readr)
  dir <- system.file("extdata", package="tximportData")
  samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
  files <- file.path(dir,"rsem", samples$run, paste0(samples$run, ".genes.results"))
  names(files) <- paste0("sample",1:6)
  txi.rsem <- tximport(files, type="rsem")
})
