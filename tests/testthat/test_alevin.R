context("alevin")
test_that("import alevin works", {

  dir <- system.file("extdata", package="tximportData")
  samps <- list.files(file.path(dir, "alevin"))
  files <- file.path(dir,"alevin",samps[1],"alevin/quants_mat.gz")
  file.exists(files)

  # regular import
  txi <- tximport(files, type="alevin")

  # no mean and variance (right now there are no boots...)
  txi.no.mv <- tximport(files, type="alevin", alevinArgs=list(dropMeanVar=TRUE))

  # again import alevin without fishpond
  txi <- tximport(files, type="alevin", alevinArgs=list(forceSlow=TRUE))

  # again import with cell barcode filtering
  txi <- tximport(files, type="alevin", alevinArgs=list(filterBarcodes=TRUE))

  # again import with tier information
  txi <- tximport(files, type="alevin", alevinArgs=list(tierImport=TRUE))
  
})
