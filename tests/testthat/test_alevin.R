context("alevin")
test_that("import alevin works", {

  dir <- system.file("extdata", package="tximportData")
  files <- file.path(dir,"alevin/neurons_900_v012/alevin/quants_mat.gz")
  file.exists(files)

  txi <- tximport(files, type="alevin")

  files <- file.path(dir,"alevin/neurons_900_v014/alevin/quants_mat.gz")
  file.exists(files)

  txi <- tximport(files, type="alevin")

  matrix.file <- file.path(dir,"alevin/neurons_900_v014/alevin/quants_mat.mtx.gz")
  mat <- Matrix::readMM(matrix.file)
  mat <- t(as.matrix(mat))
  expect_true(max(mat - unname(txi$counts)) < 1e-6)
  
})
