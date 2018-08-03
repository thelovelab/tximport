context("kallisto")
test_that("import kallisto works", {

  library(readr)
  dir <- system.file("extdata", package="tximportData")
  samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
  files <- file.path(dir,"kallisto", samples$run, "abundance.tsv.gz")
  names(files) <- paste0("sample",1:6)

  tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))

  txi <- tximport(files, type="kallisto", tx2gene=tx2gene, ignoreAfterBar=TRUE)
  expect_true(ncol(txi$counts) == length(files))
  
})
