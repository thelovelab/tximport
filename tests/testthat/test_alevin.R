context("alevin")
test_that("import alevin works", {

  dir <- system.file("extdata", package="tximportData")
  files <- file.path(dir,"alevin/neurons_900_v012/alevin/quants_mat.gz")
  file.exists(files)

  txi <- tximport(files, type="alevin")

  files <- file.path(dir,"alevin/neurons_900_v014/alevin/quants_mat.gz")
  file.exists(files)

  #txi <- tximport(files, type="alevin")
  #n <- 100
  #infrep.var <- apply(abind::abind(lapply(txi$infReps, function(x) as.matrix(x[1:n,1:n])), along=3), 1:2, var)
  #alevin.var <- as.matrix(txi$variance[1:n,1:n])
  #all.equal(alevin.var * (20/19), infrep.var, tolerance=1e-6)
  #n <- 200
  #infrep.mu <- apply(abind::abind(lapply(txi$infReps, function(x) as.matrix(x[1:n,1:n])), along=3), 1:2, mean)
  #plot(txi$counts[1:n,1:n], infrep.mu)
  
  txi <- tximport(files, type="alevin", dropInfReps=TRUE)
  idx <- 1:1000 # Bioc Windows machine can't handle the entire matrix
  cts <- unname(as.matrix(txi$counts[idx,]))

  # compare to MM import
  matrix.file <- file.path(dir,"alevin/neurons_900_v014/alevin/quants_mat.mtx.gz")
  mat <- Matrix::readMM(matrix.file)
  mat <- t(as.matrix(mat[,idx]))
  expect_true(max(abs(mat - unname(cts))) < 1e-6)

  # again import alevin without fishpond
  txi <- tximport(files, type="alevin", dropInfReps=TRUE, forceSlow=TRUE)
  idx <- 1:1000 
  cts <- unname(as.matrix(txi$counts[idx,]))
  expect_true(max(abs(mat - unname(cts))) < 1e-6)
  
})
