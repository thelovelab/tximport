context("inf reps")
test_that("inferential replicate code works", {

  library(readr)
  dir <- system.file("extdata", package="tximportData")
  samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
  files <- file.path(dir,"salmon_gibbs", samples$run, "quant.sf.gz")
  names(files) <- paste0("sample",1:6)
  
  txi <- tximport(files, type="salmon", txOut=TRUE)
  expect_true("infReps" %in% names(txi))

  txi <- tximport(files, type="salmon", txOut=TRUE, varReduce=TRUE)
  expect_true("variance" %in% names(txi))

  txi <- tximport(files, type="salmon", txOut=TRUE, dropInfReps=TRUE)
  expect_true(!any(c("infReps","variance") %in% names(txi)))
  
  # test inf replicates w/ summarization
  # (15098 txps are missing from GTF, this is Ensembl's fault, not tximport's)
  tx2gene <- read_csv(file.path(dir, "tx2gene.ensembl.v87.csv"))
  
  txi <- tximport(files, type="salmon", tx2gene=tx2gene, ignoreTxVersion=TRUE)
  expect_true(grepl("ENSG", rownames(txi$infReps[[1]])[1]))

  txi <- tximport(files, type="salmon", tx2gene=tx2gene, varReduce=TRUE, ignoreTxVersion=TRUE)
  expect_true("variance" %in% names(txi))

  txi <- tximport(files, type="salmon", tx2gene=tx2gene, dropInfReps=TRUE, ignoreTxVersion=TRUE)
  expect_true(!any(c("infReps","variance") %in% names(txi)))

  # test re-computing counts and abundances from inf replicates
  library(matrixStats)
  txi <- tximport(files, type="salmon", txOut=TRUE, infRepStat=rowMedians)
  txp <- which(rownames(txi$counts) == "ENST00000628356.2")
  expect_equal(txi$counts[txp,1], median(txi$infReps[[1]][txp,]))
  
})
