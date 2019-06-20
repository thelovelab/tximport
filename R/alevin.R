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
  con <- gzcon(file(matrix.file, "rb"))
  for (j in seq_len(num.cells)) {
    mat[,j] <- readBin(con, double(), endian = "little", n=num.genes)
  }
  close(con)
  # if inferential replicate variance exists:
  if (file.exists(var.file)) {
    counts.mat <- mat
    var.mat <- mat
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
  
  mat <- matrix(0, nrow=num.genes, ncol=num.cells, dimnames=list(gene.names, cell.names))
  con <- gzcon(file(matrix.file, "rb"))
  
  # Salmon v0.14 specific support
  num.bitvecs <- ceiling(num.genes/8)
  for (j in seq_len(num.cells)) {
    # read the bit vector
    bit.vec <- readBin(con, integer(), size=1, signed=FALSE, endian="little", n=num.bitvecs)
    bits <- sapply(bit.vec, intToBits)
    # 8 to 1, because intToBits gives the least sig bit first
    bit.ints <- apply(bits[8:1,], 2, as.integer)
    num.exp.genes <- sum(bit.ints == 1)
    # read in the expression of expressed genes
    counts <- readBin(con, double(), size=4, endian="little", n=num.exp.genes)
    locs <- head(as.vector(bit.ints), num.genes)
    mat[locs == 1,j] <- counts
  }
  close(con)
  
  mat
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
