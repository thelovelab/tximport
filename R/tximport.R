#' Import transcript-level abundances and estimated counts for gene-level analysis packages
#' 
#' @param files a character vector of filenames for the transcript-level abundances
#' @param type the type of software used to generate the abundances, 
#' which will be used to autofill the arguments below
#' @param txIn logical, whether the incoming files are transcript level (default TRUE)
#' @param txOut logical, whether the function should just output transcript-level (default FALSE)
#' @param countsFromAbundance logical, whether to generate estimated counts using 
#' abundance estimates and the average of average transcript length over samples (default FALSE). 
#' If this argument is used, then the counts are no longer correlated with average transcript length,
#' and so the length offset matrix should not be used.
#' @param geneIdCol name of column with gene id. if missing, the gene2tx argument can be used
#' @param txIdCol name of column with tx id
#' @param abundanceCol name of column with abundances (e.g. TPM or FPKM)
#' @param countsCol name of column with estimated counts
#' @param lengthCol name of column with feature length information
#' @param gene2tx a two-column data.frame linking gene id and transcript id. 
#' this is needed for software which does not provide such information in the file
#' (kallisto and Salmon)
#' @param importer a function used to read in the files
#' @param collatedFiles a character vector of filenames for software which provides
#' abundances and counts in matrix form (e.g. Cufflinks). The files should be, in order,
#' abundances, counts, and a third file with length information
#' 
#' @return a simple list with matrices: abundances, counts, length. 
#' The length matrix contains the average transcript length for each
#' gene which can be used as an offset for gene-level analysis.
#' 
#' @export
tximport <- function(files,
                     type=c("kallisto","salmon","rsem","cufflinks"),
                     txIn=TRUE,
                     txOut=FALSE,
                     countsFromAbundance=FALSE,
                     geneIdCol,
                     txIdCol,
                     abundanceCol,
                     countsCol,
                     lengthCol,
                     gene2tx=NULL,
                     importer=function(x) read.table(x,header=TRUE),
                     collatedFiles,
                     ...) {

  type <- match.arg(type, c("kallisto","salmon","rsem","cufflinks"))
  
  # kallisto presets
  if (type == "kallisto") {
    geneIdCol="gene_id"
    txIdCol <- "target_id"
    abundanceCol <- "tpm"
    countsCol <- "est_counts"
    lengthCol <- "eff_length"
    }
  
  # salmon presets
  if (type == "salmon") {
    geneIdCol="gene_id"
    txIdCol <- "Name"
    abundanceCol <- "TPM"
    countsCol <- "NumReads"
    lengthCol <- "Length"
    # because the comment lines have the same comment character as the header...
    # need to name the column names
    importer <- function(x) {
      tmp <- read.table(x,comment.char="#")
      names(tmp) <- c("Name","Length","TPM","NumReads")
      tmp
    }
  }
  
  # rsem presets
  if (type == "rsem") {
    txIn <- FALSE
    geneIdCol <- "gene_id"
    abundanceCol <- "FPKM"
    countsCol <- "expected_count"
    lengthCol <- "effective_length"
  }
  
  # if input is tx-level, need to summarize abundances, counts and lengths to gene-level
  if (txIn) {
    message("reading in files")
    for (i in seq_along(files)) {
      message(i," ",appendLF=FALSE)
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
    message("")
    
    # if the user requested just the transcript-level data:
    if (txOut) {
      return(list(abundance=abundanceMatTx, counts=countsMatTx, length=lengthMatTx))
    }
    
    # need to associate tx to genes
    # potentially remove unassociated transcript rows and warn user
    if (!is.null(gene2tx)) {
      ntxmissing <- sum(!txId %in% gene2tx$TXNAME)
      if (ntxmissing > 0) message("transcripts missing genes: ", ntxmissing)
      sub.idx <- txId %in% gene2tx$TXNAME
      abundanceMatTx <- abundanceMatTx[sub.idx,]
      countsMatTx <- countsMatTx[sub.idx,]
      lengthMatTx <- lengthMatTx[sub.idx,]
      txId <- txId[sub.idx]
      geneId <- gene2tx$GENEID[match(txId, gene2tx$TXNAME)]
    }
    
    # summarize abundance and counts
    message("summarizing abundance")
    abundanceMat <- fastby(abundanceMatTx, geneId, colSums)
    message("summarizing counts")
    countsMat <- fastby(countsMatTx, geneId, colSums)
    message("summarizing length")
    
    # the next lines calculate a weighted average of transcript length, 
    # weighting by transcript abundance.
    # this can be used as an offset / normalization factor which removes length bias
    # for the differential analysis of estimated counts summarized at the gene level.
    weightedLength <- fastby(abundanceMatTx * lengthMatTx, geneId, colSums)
    lengthMat <- weightedLength / abundanceMat   
     
    # check for NaN and if possible replace these values with geometric mean of other samples.
    # (the geometic mean here implies an offset of 0 on the log scale)
    # NaN come from samples which have abundance of 0 for all isoforms of a gene, and 
    # so we cannot calculate the weighted average. our best guess is to use the average
    # transcript length from the other samples.
    nanRows <- which(apply(lengthMat, 1, function(row) any(is.nan(row))))
    if (length(nanRows) > 0) {
      for (i in nanRows) {
        if (all(is.nan(lengthMat[i,]))) {
          lengthMat[i,] <- NA # TODO just use simple average here
        } else {
          idx <- is.nan(lengthMat[i,])
          lengthMat[i,idx] <-  exp(mean(log(lengthMat[i,!idx]), na.rm=TRUE))
        }
      }
    }
    
    if (countsFromAbundance) {
      countsSum <- colSums(countsMat)
      newCounts0 <- abundanceMat * rowMeans(lengthMat)
      newSum <- colSums(newCounts0, na.rm=TRUE)
      newCounts <- t(t(newCounts0) * (countsSum/newSum))
      countsMat <- newCounts
    }
    
    return(list(abundance=abundanceMat, counts=countsMat, length=lengthMat))
    
  # e.g. RSEM already has gene-level summaries
  # just combine the gene-level summaries across files
  } else {
    # stating the obvious:
    if (txOut) stop("txOut only an option when transcript-level data is read in (txIn=TRUE)")
  
    message("reading in files")
    for (i in seq_along(files)) {
      message(i," ",appendLF=FALSE)
      raw <- as.data.frame(importer(files[i], ...))
      stopifnot(all(c(geneIdCol, abundanceCol, lengthCol) %in% names(raw)))
      if (i == 1) {
        mat <- matrix(nrow=nrow(raw),ncol=length(files))
        rownames(mat) <- raw[[geneIdCol]]
        abundanceMat <- mat
        countsMat <- mat
        lengthMat <- mat
      }
      abundanceMat[,i] <- raw[[abundanceCol]]
      countsMat[,i] <- raw[[countsCol]]
      lengthMat[,i] <- raw[[lengthCol]]
    }
  } 
  message("")
  return(list(abundance=abundanceMat, counts=countsMat, length=lengthMat))
}

# this is much faster than by(), a bit slower than dplyr summarize_each()
fastby <- function(m, f, fun) {
  idx <- split(1:nrow(m), f)
  t(sapply(idx, function(i) fun(m[i,,drop=FALSE])))
}
