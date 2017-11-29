# summarizeToGene() splits out the summarization functions
# in tximport(), so it can be called by users to summarize
# transcript-level lists of matrices

#' @describeIn tximport Summarize transcript-level to gene-level
#' @export
summarizeToGene <- function(txi,
                            tx2gene,
                            ignoreTxVersion=FALSE,
                            countsFromAbundance=c("no","scaledTPM","lengthScaledTPM")
                            ) {

  countsFromAbundance <- match.arg(countsFromAbundance, c("no","scaledTPM","lengthScaledTPM"))
  
  # unpack matrices from list for cleaner code
  abundanceMatTx <- txi$abundance
  countsMatTx <- txi$counts
  lengthMatTx <- txi$length
  
  txId <- rownames(abundanceMatTx)
  stopifnot(all(txId == rownames(countsMatTx)))
  stopifnot(all(txId == rownames(lengthMatTx)))
  
  # need to associate tx to genes
  # potentially remove unassociated transcript rows and warn user
  if (!is.null(tx2gene)) {
    
    # code to strip dots from the rownames of matrices
    if (ignoreTxVersion) {
      txId <- sapply(strsplit(as.character(txId), "\\."), "[[", 1)
    }
    
    colnames(tx2gene) <- c("tx","gene")

    if (any(duplicated(tx2gene$tx))) {
      message("removing duplicated transcript rows from tx2gene")
      tx2gene <- tx2gene[!duplicated(tx2gene$tx),]
    }
    
    tx2gene$gene <- factor(tx2gene$gene)
    tx2gene$tx <- factor(tx2gene$tx)
    
    # if none of the rownames of the matrices (txId) are
    # in the tx2gene table something is wrong
    if (!any(txId %in% tx2gene$tx)) {
      stop("
  None of the transcripts in the quantification files are present
  in the first column of tx2gene. Check to see that you are using
  the same annotation for both.\n\n")
    }

    # remove transcripts (and genes) not in the rownames of matrices
    tx2gene <- tx2gene[tx2gene$tx %in% txId,]
    tx2gene$gene <- droplevels(tx2gene$gene)
    ntxmissing <- sum(!txId %in% tx2gene$tx)
    if (ntxmissing > 0) message("transcripts missing from tx2gene: ", ntxmissing)

    # subset to transcripts in the tx2gene table
    sub.idx <- txId %in% tx2gene$tx
    abundanceMatTx <- abundanceMatTx[sub.idx,,drop=FALSE]
    countsMatTx <- countsMatTx[sub.idx,,drop=FALSE]
    lengthMatTx <- lengthMatTx[sub.idx,,drop=FALSE]
    txId <- txId[sub.idx]

    # now create a vector of geneId which aligns to the matrices
    geneId <- tx2gene$gene[match(txId, tx2gene$tx)]
  }
  
  # summarize abundance and counts
  message("summarizing abundance")
  abundanceMat <- rowsum(abundanceMatTx, geneId)
  message("summarizing counts")
  countsMat <- rowsum(countsMatTx, geneId)
  message("summarizing length")
  
  # the next lines calculate a weighted average of transcript length, 
  # weighting by transcript abundance.
  # this can be used as an offset / normalization factor which removes length bias
  # for the differential analysis of estimated counts summarized at the gene level.
  weightedLength <- rowsum(abundanceMatTx * lengthMatTx, geneId)
  lengthMat <- weightedLength / abundanceMat   

  # pre-calculate a simple average transcript length
  # for the case the abundances are all zero for all samples.
  # first, average the tx lengths over samples
  aveLengthSamp <- rowMeans(lengthMatTx)
  # then simple average of lengths within genes (not weighted by abundance)
  aveLengthSampGene <- tapply(aveLengthSamp, geneId, mean)

  stopifnot(all(names(aveLengthSampGene) == rownames(lengthMat)))
  
  # check for NaN and if possible replace these values with geometric mean of other samples.
  # (the geometic mean here implies an offset of 0 on the log scale)
  # NaN come from samples which have abundance of 0 for all isoforms of a gene, and 
  # so we cannot calculate the weighted average. our best guess is to use the average
  # transcript length from the other samples.
  lengthMat <- replaceMissingLength(lengthMat, aveLengthSampGene)

  if (countsFromAbundance != "no") {
    countsMat <- makeCountsFromAbundance(countsMat,
                                         abundanceMat,
                                         lengthMat,
                                         countsFromAbundance)
  }
    
  return(list(abundance=abundanceMat, counts=countsMat, length=lengthMat,
              countsFromAbundance=countsFromAbundance))
}

