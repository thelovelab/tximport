context("h5")
test_that("kallisto HDF5 import works", {

  skip_if_not_installed("rhdf5")
  library(readr)
  dir <- system.file("extdata", package="tximportData")
  samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
  tsv <- file.path(dir,"kallisto", samples$run, "abundance.tsv.gz")

  # tximportData no longer ships kallisto abundance.h5 files, so write
  # small ones (with bootstraps) from the kallisto TSV output
  n <- 1000
  nboot <- 5
  files <- file.path(tempdir(), "kallisto_h5", samples$run, "abundance.h5")
  for (i in seq_along(files)) {
    dir.create(dirname(files[i]), recursive=TRUE, showWarnings=FALSE)
    unlink(files[i])
    ab <- read_tsv(tsv[i], n_max=n, show_col_types=FALSE)
    rhdf5::h5createFile(files[i])
    rhdf5::h5createGroup(files[i], "aux")
    rhdf5::h5createGroup(files[i], "bootstrap")
    rhdf5::h5write(ab$target_id, files[i], "aux/ids")
    rhdf5::h5write(ab$length, files[i], "aux/lengths")
    rhdf5::h5write(ab$eff_length, files[i], "aux/eff_lengths")
    rhdf5::h5write(ab$est_counts, files[i], "est_counts")
    for (b in seq_len(nboot) - 1) {
      rhdf5::h5write(ab$est_counts * runif(n, .9, 1.1), files[i],
                     paste0("bootstrap/bs", b))
    }
  }
  names(files) <- paste0("sample",1:2)

  txi <- tximport(files, type="kallisto", txOut=TRUE)
  expect_true("infReps" %in% names(txi))
  expect_equal(dim(txi$infReps$sample1), c(n, nboot))
  txi <- tximport(files, type="kallisto", txOut=TRUE, varReduce=TRUE)
  expect_true("variance" %in% names(txi))
  txi <- tximport(files, type="kallisto", txOut=TRUE, dropInfReps=TRUE)
  expect_true(!any(c("infReps","variance") %in% names(txi)))

})
