readAlevinPreV014 <- function(files) {
  message("using importer for pre-v0.14.0 Alevin files")
  dir <- sub("/alevin$","",dirname(files))
  barcode.file <- file.path(dir, "alevin/quants_mat_rows.txt")
  gene.file <- file.path(dir, "alevin/quants_mat_cols.txt")
  matrix.file <- file.path(dir, "alevin/quants_mat.gz")
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
    message("reading in alevin inferential variance")
    con <- gzcon(file(var.file, "rb"))
    for (j in seq_len(num.cells)) {
      var.mat[,j] <- readBin(con, double(), endian = "little", n=num.genes)
    }
    close(con)
    mat <- list(counts.mat, var.mat)
  }
  mat
}

readAlevin <- function(files) {
  dir <- sub("/alevin$","",dirname(files))  
  barcode.file <- file.path(dir, "alevin/quants_mat_rows.txt")
  gene.file <- file.path(dir, "alevin/quants_mat_cols.txt")
  matrix.file <- file.path(dir, "alevin/quants_mat.gz")
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

  message("reading in alevin gene-level counts across cells")
  mat <- readAlevinBits(matrix.file, gene.names, cell.names)
  
  if (num.boot > 0) {

    message("reading in alevin inferential variance")
    var.mat <- readAlevinBits(var.file, gene.names, cell.names)
    return(list(mat, var.mat))
    
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

readAlevinBits <- function(matrix.file, gene.names, cell.names) {
  num.cells <- length(cell.names)
  num.genes <- length(gene.names)
  len.bit.vec <- ceiling(num.genes/8)
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
                              dimnames=list(gene.names, cell.names))
  mat
}
