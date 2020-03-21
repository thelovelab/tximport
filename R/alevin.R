readAlevinPreV014 <- function(files) {
  message("using importer for pre-v0.14.0 Alevin files")
  message("pre-v0.14.0 Alevin import will be deprecated in October 2020")
  dir <- sub("/alevin$","",dirname(files))
  barcode.file <- file.path(dir, "alevin/quants_mat_rows.txt")
  gene.file <- file.path(dir, "alevin/quants_mat_cols.txt")
  matrix.file <- file.path(dir, "alevin/quants_mat.gz")
  boot.barcode.file <- file.path(dir, "alevin/quants_boot_rows.txt")
  var.file <- file.path(dir, "alevin/quants_var_mat.gz")
  for (f in c(barcode.file, gene.file, matrix.file)) {
    if (!file.exists(f)) {
      stop("expecting 'files' to point to 'quants_mat.gz' file in a directory 'alevin'
  also containing 'quants_mat_rows.txt' and 'quant_mat_cols.txt'.
  please re-run alevin preserving output structure")
    }
  }
  cell.names <- readLines(barcode.file)
  gene.names <- readLines(gene.file)
  num.cells <- length(cell.names)
  num.genes <- length(gene.names)
  mat <- matrix(nrow=num.genes, ncol=num.cells, dimnames=list(gene.names, cell.names))
  message("reading in alevin gene-level counts across cells")
  con <- gzcon(file(matrix.file, "rb"))
  for (j in seq_len(num.cells)) {
    mat[,j] <- readBin(con, double(), endian = "little", n=num.genes)
  }
  close(con)
  # if inferential replicate variance exists:
  if (file.exists(var.file)) {
    counts.mat <- mat
    var.mat <- mat
    # need to re-arrange to match the counts matrix
    suppressWarnings({
      boot.cell.names <- readLines(boot.barcode.file)[-1]
    })
    stopifnot(length(boot.cell.names) == length(cell.names))
    colnames(var.mat) <- boot.cell.names
    message("reading in alevin inferential variance")
    con <- gzcon(file(var.file, "rb"))
    for (j in seq_len(num.cells)) {
      var.mat[,j] <- readBin(con, double(), endian = "little", n=num.genes)
    }
    var.mat <- var.mat[,cell.names]
    close(con)
    mat <- list(counts=counts.mat, variance=var.mat)
  }
  mat
}

readAlevin <- function(files, dropInfReps, forceSlow) {
  dir <- sub("/alevin$","",dirname(files))  
  barcode.file <- file.path(dir, "alevin/quants_mat_rows.txt")
  gene.file <- file.path(dir, "alevin/quants_mat_cols.txt")
  matrix.file <- file.path(dir, "alevin/quants_mat.gz")
  mean.file <- file.path(dir, "alevin/quants_mean_mat.gz")
  var.file <- file.path(dir, "alevin/quants_var_mat.gz")
  boot.file <- file.path(dir, "alevin/quants_boot_mat.gz")
  boot.barcode.file <- file.path(dir, "alevin/quants_boot_rows.txt")
  for (f in c(barcode.file, gene.file, matrix.file)) {
    if (!file.exists(f)) {
      stop("expecting 'files' to point to 'quants_mat.gz' file in a directory 'alevin'
  also containing 'quants_mat_rows.txt' and 'quant_mat_cols.txt'.
  please re-run alevin preserving output structure")
    }
  }
  cell.names <- readLines(barcode.file)
  gene.names <- readLines(gene.file)
  
  if (!requireNamespace("jsonlite", quietly=TRUE)) {
    stop("importing alevin requires package `jsonlite`")
  }
  jsonPath <- file.path(dir, "cmd_info.json")
  cmd_info <- jsonlite::fromJSON(jsonPath)
  if ("numCellBootstraps" %in% names(cmd_info)) {
    num.boot <- as.numeric(cmd_info$numCellBootstraps)
  } else {
    num.boot <- 0
  }

  if (!requireNamespace("Matrix", quietly=TRUE)) {
    stop("importing alevin requires package `Matrix`")
  }

  # test for fishpond >= 1.1.17
  hasFishpond <- TRUE
  if (!requireNamespace("fishpond", quietly=TRUE)) {
    hasFishpond <- FALSE
  } else {
    if (packageVersion("fishpond") < "1.1.18") {
      hasFishpond <- FALSE
    }
  }
  if (!hasFishpond)
    message("importing alevin data is much faster after installing `fishpond` (>= 1.1.18)")

  # for testing purposes, force the use of the slower R code for importing alevin
  if (forceSlow) {
    hasFishpond <- FALSE
  }
  
  extraMsg <- if (hasFishpond) "with fishpond" else ""
  message(paste("reading in alevin gene-level counts across cells", extraMsg))

  if (hasFishpond) {
    # reads alevin's Efficient Data Storage (EDS) format
    # using C++ code in the fishpond package
    mat <- readAlevinFast(matrix.file, gene.names, cell.names)
  } else {
    # reads alevin EDS format in R, using e.g. `readBin` and `intToBits`
    # slow in R, because requires looping over cells to read positions and expression
    mat <- readAlevinBits(matrix.file, gene.names, cell.names)
  }
  
  if (num.boot > 0) {

    message(paste("reading in alevin inferential variance", extraMsg))
    mean.exists <- file.exists(mean.file)
    var.exists <- file.exists(var.file)
    boot.exists <- file.exists(boot.file)
    stopifnot(mean.exists)
    stopifnot(var.exists)

    boot.cell.names <- readLines(boot.barcode.file)

    if (hasFishpond) {
      mean.mat <- readAlevinFast(mean.file, gene.names, boot.cell.names)
      var.mat <- readAlevinFast(var.file, gene.names, boot.cell.names)
    } else {
      mean.mat <- readAlevinBits(mean.file, gene.names, boot.cell.names)
      var.mat <- readAlevinBits(var.file, gene.names, boot.cell.names)
    }

    # need to re-arrange to match the counts matrix
    mean.mat <- mean.mat[,cell.names]
    var.mat <- var.mat[,cell.names]
    
    if (boot.exists & !dropInfReps) {
      # read in bootstrap inferential replicates
      message("reading in alevin inferential replicates (set dropInfReps=TRUE to skip)")
      infReps <- readAlevinInfReps(boot.file, gene.names, boot.cell.names, num.boot)
      # need to re-arrange to match the counts matrix
      infReps <- lapply(infReps, function(z) z[,cell.names])
      return(list(counts=mat, mean=mean.mat, variance=var.mat, infReps=infReps))
    } else {
      # return counts, mean and variance
      return(list(counts=mat, mean=mean.mat, variance=var.mat))
    }
    
  } else {
    
    return(mat)
    
  }
}

getAlevinVersion <- function(files) {
  if (!requireNamespace("jsonlite", quietly=TRUE)) {
    stop("importing Alevin quantification requires package `jsonlite`")
  }
  fish_dir <- dirname(dirname(files))
  jsonPath <- file.path(fish_dir, "cmd_info.json")
  cmd_info <- jsonlite::fromJSON(jsonPath)
  cmd_info$salmon_version
}

# this is the R version of the reader for alevin's EDS format,
# see below for another function that leverages C++ code from fishpond::readEDS()
readAlevinBits <- function(matrix.file, gene.names, cell.names) {
  num.cells <- length(cell.names)
  num.genes <- length(gene.names)
  len.bit.vec <- ceiling(num.genes/8)
  # the bit vector matrix is 8 rows x num.genes/8 columns
  # and stores the positions of the non-zero counts
  bits.mat <- matrix(nrow=8, ncol=num.cells * len.bit.vec)
  con <- gzcon(file(matrix.file, "rb"))
  for (j in seq_len(num.cells)) {
    # read the bit vectors
    ints <- readBin(con, integer(), size=1, signed=FALSE, endian="little", n=len.bit.vec)
    bits <- matrix(intToBits(ints), nrow=32)
    mode(bits) <- "integer"
    # 8 to 1, because intToBits gives the least sig bit first
    bits <- bits[8:1,]
    num.exp.genes <- sum(bits == 1)
    # store bits in the matrix
    idx <- (j-1) * len.bit.vec + seq_len(len.bit.vec)
    bits.mat[,idx] <- bits
    # read in counts, but don't store
    counts <- readBin(con, double(), size=4, endian="little", n=num.exp.genes)
  }
  close(con)
  con <- gzcon(file(matrix.file, "rb"))

  # stores all of the non-zero counts from all cells concatenated
  counts.vec <- numeric(sum(bits.mat))
  ptr <- 0
  for (j in seq_len(num.cells)) {
    # read in bit vectors, but don't store
    ints <- readBin(con, integer(), size=1, signed=FALSE, endian="little", n=len.bit.vec)
    bit.idx <- (j-1) * len.bit.vec + seq_len(len.bit.vec)
    num.exp.genes <- sum(bits.mat[,bit.idx])
    cts.idx <- ptr + seq_len(num.exp.genes)
    counts.vec[cts.idx] <- readBin(con, double(), size=4, endian="little", n=num.exp.genes)
    ptr <- ptr + num.exp.genes
  }
  close(con)

  gene.idx <- lapply(seq_len(num.cells), function(j) {
    idx <- (j-1) * len.bit.vec + seq_len(len.bit.vec)
    which(head(as.vector(bits.mat[,idx]), num.genes) == 1)
  })
  len.gene.idx <- lengths(gene.idx)
  cell.idx <- rep(seq_along(len.gene.idx), len.gene.idx)

  # build sparse matrix
  mat <- Matrix::sparseMatrix(i=unlist(gene.idx),
                              j=cell.idx,
                              x=counts.vec,
                              dims=c(num.genes, num.cells),
                              dimnames=list(gene.names, cell.names),
                              giveCsparse=FALSE)
  mat
}

# this function performs the same operation as the above R code,
# reading in alevin's EDS format and creating a sparse matrix,
# but it leverages the C++ code in fishpond::readEDS()
readAlevinFast <- function(matrix.file, gene.names, cell.names) {
  num.cells <- length(cell.names)
  num.genes <- length(gene.names)
  mat <- fishpond::readEDS(num.genes, num.cells, matrix.file)
  dimnames(mat) <- list(gene.names, cell.names)
  mat
}

readAlevinInfReps <- function(boot.file, gene.names, cell.names, num.boot) {
  num.cells <- length(cell.names)
  num.genes <- length(gene.names)
  len.bit.vec <- ceiling(num.genes/8)
  # a list to store the bit vector matrices
  # each element of the list is for an inf rep
  # (for description of `bits.mat` see the readAlevinBits() function)
  bits.mat.list <- lapply(seq_len(num.boot), function(i)
    matrix(nrow=8, ncol=num.cells * len.bit.vec))
  con <- gzcon(file(boot.file, "rb"))
  for (j in seq_len(num.cells)) {
    for (i in seq_len(num.boot)) {
      # read the bit vectors
      ints <- readBin(con, integer(), size=1, signed=FALSE, endian="little", n=len.bit.vec)
      bits <- matrix(intToBits(ints), nrow=32)
      mode(bits) <- "integer"
      # 8 to 1, because intToBits gives the least sig bit first
      bits <- bits[8:1,]
      num.exp.genes <- sum(bits == 1)
      # store bits in the matrix
      idx <- (j-1) * len.bit.vec + seq_len(len.bit.vec)
      bits.mat.list[[i]][,idx] <- bits
      # read in counts, but don't store
      counts <- readBin(con, double(), size=4, endian="little", n=num.exp.genes)
    }
  }
  close(con)
  con <- gzcon(file(boot.file, "rb"))

  # as with `bits.mat.list`, `counts.vec.list` is a list
  # over inferential replicates, each one is a counts vector
  # storing all of the non-zero counts from all cells concatenated
  counts.vec.list <- lapply(seq_len(num.boot), function(i)
    numeric(sum(bits.mat.list[[i]])))
  ptr <- numeric(num.boot)
  for (j in seq_len(num.cells)) {
    for (i in seq_len(num.boot)) {
      # read in bit vectors, but don't store
      ints <- readBin(con, integer(), size=1, signed=FALSE, endian="little", n=len.bit.vec)
      bit.idx <- (j-1) * len.bit.vec + seq_len(len.bit.vec)
      num.exp.genes <- sum(bits.mat.list[[i]][,bit.idx])
      cts.idx <- ptr[i] + seq_len(num.exp.genes)
      counts.vec.list[[i]][cts.idx] <- readBin(con, double(), size=4, endian="little", n=num.exp.genes)
      ptr[i] <- ptr[i] + num.exp.genes
    }
  }
  close(con)

  infReps <- lapply(seq_len(num.boot), function(i) {
    gene.idx <- lapply(seq_len(num.cells), function(j) {
      idx <- (j-1) * len.bit.vec + seq_len(len.bit.vec)
      which(head(as.vector(bits.mat.list[[i]][,idx]), num.genes) == 1)
    })
    len.gene.idx <- lengths(gene.idx)
    cell.idx <- rep(seq_along(len.gene.idx), len.gene.idx)
    z <- Matrix::sparseMatrix(i=unlist(gene.idx),
                         j=cell.idx,
                         x=counts.vec.list[[i]],
                         dims=c(num.genes, num.cells),
                         dimnames=list(gene.names, cell.names),
                         giveCsparse=FALSE)
  })

  infReps
}
