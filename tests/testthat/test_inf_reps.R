context("inf reps")
test_that("inferential replicate code works", {

  library(readr)
  dir <- system.file("extdata", package="tximportData")
  samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)

  # tximportData no longer ships Gibbs samples, so copy the salmon
  # output and write synthetic Gibbs samples in the salmon format
  ngibbs <- 5
  files <- file.path(tempdir(), "salmon_gibbs", samples$run, "quant.sf.gz")
  for (i in seq_along(files)) {
    out <- dirname(files[i])
    dir.create(file.path(out, "aux_info", "bootstrap"),
               recursive=TRUE, showWarnings=FALSE)
    file.copy(file.path(dir, "salmon", samples$run[i], "quant.sf.gz"), out,
              overwrite=TRUE)
    file.copy(file.path(dir, "salmon", samples$run[i], "cmd_info.json"), out,
              overwrite=TRUE)
    minfo <- jsonlite::fromJSON(file.path(dir, "salmon", samples$run[i],
                                          "aux_info", "meta_info.json"))
    minfo$samp_type <- "gibbs"
    minfo$num_bootstraps <- ngibbs
    jsonlite::write_json(minfo, file.path(out, "aux_info", "meta_info.json"),
                         auto_unbox=TRUE, pretty=TRUE)
    quant <- read_tsv(files[i], show_col_types=FALSE)
    stopifnot(nrow(quant) == minfo$num_targets)
    gibbs <- quant$NumReads * runif(nrow(quant) * ngibbs, .9, 1.1)
    con <- gzfile(file.path(out, "aux_info", "bootstrap", "bootstraps.gz"), "wb")
    writeBin(gibbs, con)
    close(con)
  }
  names(files) <- paste0("sample",1:2)

  txi <- tximport(files, type="salmon", txOut=TRUE)
  expect_true("infReps" %in% names(txi))
  expect_equal(dim(txi$infReps$sample1), c(nrow(txi$counts), ngibbs))

  txi <- tximport(files, type="salmon", txOut=TRUE, varReduce=TRUE)
  expect_true("variance" %in% names(txi))

  txi <- tximport(files, type="salmon", txOut=TRUE, dropInfReps=TRUE)
  expect_true(!any(c("infReps","variance") %in% names(txi)))

  # test inf replicates w/ summarization
  tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))

  txi <- tximport(files, type="salmon", tx2gene=tx2gene)
  expect_true(grepl("ENSG", rownames(txi$infReps[[1]])[1]))

  txi <- tximport(files, type="salmon", tx2gene=tx2gene, varReduce=TRUE)
  expect_true("variance" %in% names(txi))

  txi <- tximport(files, type="salmon", tx2gene=tx2gene, dropInfReps=TRUE)
  expect_true(!any(c("infReps","variance") %in% names(txi)))

  # test re-computing counts and abundances from inf replicates
  library(matrixStats)
  txi <- tximport(files, type="salmon", txOut=TRUE, infRepStat=rowMedians)
  txp <- which.max(txi$counts[,1])
  expect_equal(txi$counts[txp,1], median(txi$infReps[[1]][txp,]))

})
