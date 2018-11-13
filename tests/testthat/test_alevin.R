context("alevin")
test_that("import alevin works", {

  dir <- system.file("extdata", package="tximportData")
  files <- file.path(dir,"alevin/neurons_900/alevin/quants_mat.gz")
  file.exists(files)

  txi <- tximport(files, type="alevin")
  names(txi)
  dim(txi$counts)

})
