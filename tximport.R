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
  
  # if tx-level, need to summarize to gene-level
  if (level == "tx") {
    cat("reading in files: ")
    for (i in seq_along(files)) {
      cat(i,"")
      raw <- as.data.frame(importer(files[i], ...))
      # does the table contain gene association or was an external gene2tx table provided?
      if (is.null(gene2tx)) {
        # Cufflinks includes the gene ID in the table
        stopifnot(all(c(geneIdCol, lengthCol, abundanceCol) %in% names(raw)))
        if (i == 1) {
          geneId <- raw[[geneIdCol]]
        } else {
          stopifnot(all(geneId == raw[[geneIdCol]]))
        }
      } else {
        # Salmon and kallisto do not include the gene ID, need an external table
        stopifnot(all(c(lengthCol, abundanceCol) %in% names(raw)))
        if (i == 1) {
          txId <- raw[[txIdCol]]
        } else {
          stopifnot(all(txId == raw[[txIdCol]]))
        }
      }
      # allocate matrices
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
    
    # need to associate tx to genes, and remove unassociated rows
    if (!is.null(gene2tx)) {
      cat("\ntranscripts missing genes:",sum(!txId %in% gene2tx$TXNAME))
      sub.idx <- txId %in% gene2tx$TXNAME
      abundanceMatTx <- abundanceMatTx[sub.idx,]
      countsMatTx <- countsMatTx[sub.idx,]
      lengthMatTx <- lengthMatTx[sub.idx,]
      txId <- txId[sub.idx]
      geneId <- gene2tx$GENEID[match(txId, gene2tx$TXNAME)]
    }
    
    cat("\nsummarizing abundance")
    abundanceMat <- fastby(abundanceMatTx, geneId, colSums)
    cat("\nsummarizing counts")
    countsMat <- fastby(countsMatTx, geneId, colSums)
    cat("\nsummarizing length\n")
    weightedLength <- fastby(abundanceMatTx * lengthMatTx, geneId, colSums)
    # a weighted average of transcript length, weighting by transcript abundance
    lengthMat <- weightedLength / abundanceMat   
 
    return(list(abundance=abundanceMat, counts=countsMat, length=lengthMat))
    
  # e.g. RSEM already has gene-level summaries
  # just cbind() the gene-level summary across files
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


  