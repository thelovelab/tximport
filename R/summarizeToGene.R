#' @rdname summarizeToGene
#' @export
setGeneric("summarizeToGene", function(object, ...) standardGeneric("summarizeToGene"))

summarizeToGene.list <- function(object,
                                 tx2gene,
                                 varReduce=FALSE,
                                 ignoreTxVersion=FALSE,
                                 ignoreAfterBar=FALSE,
                                 countsFromAbundance=c("no","scaledTPM","lengthScaledTPM")
                                 ) {

  countsFromAbundance <- match.arg(countsFromAbundance, c("no","scaledTPM","lengthScaledTPM"))

  if (!is.null(object$countsFromAbundance)) {
    if (countsFromAbundance == "no" & object$countsFromAbundance != "no") {
      warning(paste0("the incoming counts have countsFromAbundance = '",
                     object$countsFromAbundance,"',
  and so the original counts are no longer accessible.
  to use countsFromAbundance='no', re-run objectmport() with this setting.
  over-riding 'countsFromAbundance' to set it to: ",
  object$countsFromAbundance))
      countsFromAbundance <- object$countsFromAbundance
    }
  }
  
  # unpack matrices from list for cleaner code
  abundanceMatTx <- object$abundance
  countsMatTx <- object$counts
  lengthMatTx <- object$length
  
  txId <- rownames(abundanceMatTx)
  stopifnot(all(txId == rownames(countsMatTx)))
  stopifnot(all(txId == rownames(lengthMatTx)))
  
  # need to associate tx to genes
  # potentially remove unassociated transcript rows and warn user
  if (!is.null(tx2gene)) {

    # code to strip dots or bars and all remaining chars from the rownames of matrices
    if (ignoreTxVersion) {
      txId <- sub("\\..*", "", txId)
    } else if (ignoreAfterBar) {
      txId <- sub("\\|.*", "", txId)
    }
    
    tx2gene <- cleanTx2Gene(tx2gene)    
    
    # if none of the rownames of the matrices (txId) are
    # in the tx2gene table something is wrong
    if (!any(txId %in% tx2gene$tx)) {
      txFromFile <- paste0("Example IDs (file): [", paste(head(txId,3),collapse=", "),", ...]")
      txFromTable <- paste0("Example IDs (tx2gene): [", paste(head(tx2gene$tx,3),collapse=", "),", ...]")
      stop(paste0("
  None of the transcripts in the quantification files are present
  in the first column of tx2gene. Check to see that you are using
  the same annotation for both.\n\n",txFromFile,"\n\n",txFromTable,
  "\n\n  This can sometimes (not always) be fixed using 'ignoreTxVersion' or 'ignoreAfterBar'.\n\n"))
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

  if ("infReps" %in% names(object)) {
    infReps <- lapply(object$infReps, function(x) rowsum(x[sub.idx,,drop=FALSE], geneId))
    message("summarizing inferential replicates")
  }
  
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


  if ("infReps" %in% names(object)) {
    if (varReduce) {
      vars <- sapply(infReps, rowVars)
      out <- list(abundance=abundanceMat,
                  counts=countsMat, variance=vars,
                  length=lengthMat,
                  countsFromAbundance=countsFromAbundance)
    } else {
      out <- list(abundance=abundanceMat,
                  counts=countsMat, infReps=infReps,
                  length=lengthMat,
                  countsFromAbundance=countsFromAbundance)
    }
  } else {
    out <- list(abundance=abundanceMat,
                counts=countsMat,
                length=lengthMat,
                countsFromAbundance=countsFromAbundance)
  }
  
  return(out)
}

#' Summarize estimated quantitites to gene-level
#'
#' Summarizes abundances, counts, lengths, (and inferential
#' replicates or variance) from transcript- to gene-level.
#'
#' @param object the list of matrices of trancript-level abundances,
#' counts, lengths produced by \code{\link{tximport}},
#' with a \code{countsFromAbundance} element that tells
#' how the counts were generated.
#' @param tx2gene see \code{\link{tximport}}
#' @param varReduce see \code{\link{tximport}}
#' @param ignoreTxVersion see \code{\link{tximport}}
#' @param ignoreAfterBar see \code{\link{tximport}}
#' @param countsFromAbundance see \code{\link{tximport}}
#' @param ... additional arguments, ignored
#'
#' @return a list of matrices of gene-level abundances, counts, lengths,
#' (and inferential replicates or variance if inferential replicates
#' are present).
#'
#' @rdname summarizeToGene
#' @docType methods
#' @aliases summarizeToGene,list-method
#'
#' @seealso \code{\link{tximport}}
#' 
#' @export
setMethod("summarizeToGene", signature(object="list"),
          summarizeToGene.list)

cleanTx2Gene <- function(tx2gene) {
  colnames(tx2gene) <- c("tx","gene")
  if (any(duplicated(tx2gene$tx))) {
    message("removing duplicated transcript rows from tx2gene")
    tx2gene <- tx2gene[!duplicated(tx2gene$tx),]
  }
  tx2gene$gene <- factor(tx2gene$gene)
  tx2gene$tx <- factor(tx2gene$tx)
  tx2gene
}
