# this is much faster than by(), a bit slower than dplyr summarize_each()
fastby <- function(m, f, fun) {
  idx <- split(1:nrow(m), f)
  t(sapply(idx, function(i) fun(m[i,,drop=FALSE])))
}

tximport <- function(files,
                     level=c("tx","gene"),
                     geneIdCol="gene_id",
                     txIdCol="target_id",
                     lengthCol="eff_length",
                     abundanceCol="tpm",
                     countsCol="est_counts",
                     gene2tx=NULL,
                     importer=function(x) read.table(x,header=TRUE),
                     ...) {

  # this is the level of the input files
  level <- match.arg(level, c("tx","gene"))
  
  # if input is tx-level, need to summarize to gene-level
  if (level == "tx") {
    cat("reading in files: ")
    for (i in seq_along(files)) {
      cat(i,"")
      raw <- as.data.frame(importer(files[i], ...))
      # does the table contain gene association or was an external gene2tx table provided?
      if (is.null(gene2tx)) {
        # e.g. Cufflinks includes the gene ID in the table
        stopifnot(all(c(geneIdCol, lengthCol, abundanceCol) %in% names(raw)))
        if (i == 1) {
          geneId <- raw[[geneIdCol]]
        } else {
          stopifnot(all(geneId == raw[[geneIdCol]]))
        }
      } else {
        # e.g. Salmon and kallisto do not include the gene ID, need an external table
        stopifnot(all(c(lengthCol, abundanceCol) %in% names(raw)))
        if (i == 1) {
          txId <- raw[[txIdCol]]
        } else {
          stopifnot(all(txId == raw[[txIdCol]]))
        }
      }
      # create empty matrices
      if (i == 1) {
        mat <- matrix(nrow=nrow(raw),ncol=length(files))
        abundanceMatTx <- mat
        countsMatTx <- mat
        lengthMatTx <- mat
      }
      abundanceMatTx[,i] <- raw[[abundanceCol]]
      countsMatTx[,i] <- raw[[countsCol]]
      lengthMatTx[,i] <- raw[[lengthCol]]
    }
    # need to associate tx to genes, and remove unassociated rows and warn user
    if (!is.null(gene2tx)) {
      cat("\ntranscripts missing genes:",sum(!txId %in% gene2tx$TXNAME))
      sub.idx <- txId %in% gene2tx$TXNAME
      abundanceMatTx <- abundanceMatTx[sub.idx,]
      countsMatTx <- countsMatTx[sub.idx,]
      lengthMatTx <- lengthMatTx[sub.idx,]
      txId <- txId[sub.idx]
      geneId <- gene2tx$GENEID[match(txId, gene2tx$TXNAME)]
    }
    # summarize abundance and counts
    cat("\nsummarizing abundance")
    abundanceMat <- fastby(abundanceMatTx, geneId, colSums)
    cat("\nsummarizing counts")
    countsMat <- fastby(countsMatTx, geneId, colSums)
    cat("\nsummarizing length\n")
    
    # the next two lines, calculate a weighted average of transcript length, 
    # weighting by transcript abundance.
    # this can be used as an offset / normalization factor which removes length bias
    # for the differential analysis of estimated counts summarized at the gene level
    weightedLength <- fastby(abundanceMatTx * lengthMatTx, geneId, colSums)
    lengthMat <- weightedLength / abundanceMat   
 
    # check for NaN and if possible replace these values with geometric mean of other samples.
    # NaN come from samples which have abundance of 0 for all isoforms of a gene, and 
    # so we cannot calculate the weighted average. our best guess is to use the average
    # transcript length from the other samples.
    nanRows <- which(apply(lengthMat, 1, function(row) any(is.nan(row))))
    if (length(nanRows) > 0) {
      for (i in nanRows) {
        if (all(is.nan(lengthMat[i,]))) {
          lengthMat[i,] <- NA
        } else {
          idx <- is.nan(lengthMat[i,])
          lengthMat[i,idx] <-  exp(mean(log(lengthMat[i,!idx]), na.rm=TRUE))
        }
      }
    }
    
    return(list(abundance=abundanceMat, counts=countsMat, length=lengthMat))
    
  # e.g. RSEM already has gene-level summaries
  # just combine the gene-level summaries across files
  } else if (level == "gene") {
    cat("reading in files: ")
    for (i in seq_along(files)) {
      cat(i,"")
      raw <- as.data.frame(importer(files[i], ...))
      stopifnot(all(c(geneIdCol, abundanceCol, lengthCol) %in% names(raw)))
      if (i == 1) {
        mat <- matrix(nrow=nrow(raw),ncol=length(files))
        abundanceMatTx <- mat
        countsMatTx <- mat
        lengthMatTx <- mat
      }
      abundanceMatTx[,i] <- raw[[abundanceCol]]
      countsMatTx[,i] <- raw[[countsCol]]
      lengthMatTx[,i] <- raw[[lengthCol]]
    }
  } 
  cat("\n")
  return(list(abundance=abundanceMat, counts=countsMat, length=lengthMat))
}

