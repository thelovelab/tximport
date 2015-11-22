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
  
  lengthList <- list()
  abundanceList <- list()
  countsList <- list()
  
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
      lengthList[[i]] <- raw[[lengthCol]]
      abundanceList[[i]] <- raw[[abundanceCol]]
      countsList[[i]] <- raw[[countsCol]]
    }
    
    lengthMatTx <- as.matrix(do.call(cbind, lengthList))
    abundanceMatTx <- as.matrix(do.call(cbind, abundanceList))
    countsMatTx <- as.matrix(do.call(cbind, countsList))
    
    # need to associate tx to genes, and remove unassociated rows
    if (!is.null(gene2tx)) {
      cat("\ntx missing genes:",sum(!txId %in% gene2tx$TXNAME))
      sub.idx <- txId %in% gene2tx$TXNAME
      lengthMatTx <- lengthMatTx[sub.idx,]
      abundanceMatTx <- abundanceMatTx[sub.idx,]
      countsMatTx <- countsMatTx[sub.idx,]
      txId <- txId[sub.idx]
      geneId <- gene2tx$GENEID[match(txId, gene2tx$TXNAME)]
    }
    
    cat("\nsummarizing abundance")
    abundanceMat <- do.call(rbind, by(abundanceMatTx, geneId, colSums))
    cat("\nsummarizing counts")
    countsMat <- do.call(rbind, by(countsMatTx, geneId, colSums))
    cat("\nsummarizing length\n")
    weightedLength <- do.call(rbind, by(abundanceMatTx * lengthMatTx, geneId, colSums))
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
      abundanceList[[i]] <- raw[[abundanceCol]]
      names(abundanceList[[i]]) <- raw[[geneIdCol]]
      countsList[[i]] <- raw[[countsCol]]
      names(countsList[[i]]) <- raw[[geneIdCol]]
      lengthList[[i]] <- raw[[lengthCol]]
      names(lengthList[[i]]) <- raw[[geneIdCol]]
    }
    abundanceMat <- as.matrix(do.call(cbind, abundanceList))
    countsMat <- as.matrix(do.call(cbind, countsList))
    lengthMat <- as.matrix(do.call(cbind, lengthList))
  } 
  cat("\n")
  return(list(abundance=abundanceMat, counts=countsMat, length=lengthMat))
}


  