#' Import transcript-level abundances and estimated counts for gene-level analysis packages
#'
#' \code{tximport} imports transcript-level estimates from various
#' external software and optionally summarizes abundances, counts, and transcript lengths
#' to the gene-level (default) or outputs transcript-level matrices
#' (see \code{txOut} argument).
#' While \code{tximport} summarizes to the gene-level by default, 
#' the user can also perform the import and summarization steps manually,
#' by specifing \code{txOut=TRUE} and then using the function \code{summarizeToGene}.
#' Note however that this is equivalent to \code{tximport} with
#' \code{txOut=FALSE} (the default).
#'
#' \strong{Solutions} to the error "tximport failed at summarizing to the gene-level":
#'
#' \enumerate{
#'   \item provide a \code{tx2gene} data.frame linking transcripts to genes (more below)
#'   \item avoid gene-level summarization by specifying \code{txOut=TRUE}
#'   \item set \code{geneIdCol} to an appropriate column in the files
#' }
#' 
#' See \code{vignette('tximport')} for example code for generating a
#' \code{tx2gene} data.frame from a \code{TxDb} object.
#' Note that the \code{keys} and \code{select} functions used
#' to create the \code{tx2gene} object are documented
#' in the man page for \link[AnnotationDbi]{AnnotationDb-class} objects
#' in the AnnotationDbi package (TxDb inherits from AnnotationDb).
#' For further details on generating TxDb objects from various inputs
#' see \code{vignette('GenomicFeatures')} from the GenomicFeatures package.
#'
#' \strong{Version support}: The last known supported versions of the
#' external quantifiers are:
#' kallisto 0.42.4, Salmon 0.6.0, Sailfish 0.9.0, RSEM 1.2.11.
#'
#' @param files a character vector of filenames for the transcript-level abundances
#' @param type character, the type of software used to generate the abundances.
#' Options are "kallisto", "salmon", "sailfish", "rsem".
#' This argument is used to autofill the arguments below (geneIdCol, etc.)
#' "none" means that the user will specify these columns.
#' @param txIn logical, whether the incoming files are transcript level (default TRUE)
#' @param txOut logical, whether the function should just output transcript-level (default FALSE)
#' @param countsFromAbundance character, either "no" (default), "scaledTPM", or "lengthScaledTPM",
#' for whether to generate estimated counts using abundance estimates scaled up to library size
#' (scaledTPM) or additionally scaled using the average transcript length over samples and
#' the library size (lengthScaledTPM). if using scaledTPM or lengthScaledTPM, 
#' then the counts are no longer correlated with average transcript length,
#' and so the length offset matrix should not be used.
#' @param tx2gene a two-column data.frame linking transcript id (column 1) to gene id (column 2).
#' the column names are not relevant, but this column order must be used. 
#' this argument is required for gene-level summarization for methods
#' that provides transcript-level estimates only
#' (kallisto, Salmon, Sailfish)
#' @param reader a function to replace read.delim in the pre-set importer functions,
#' for example substituting read_tsv from the readr package will substantially 
#' speed up tximport
#' @param geneIdCol name of column with gene id. if missing, the gene2tx argument can be used
#' @param txIdCol name of column with tx id
#' @param abundanceCol name of column with abundances (e.g. TPM or FPKM)
#' @param countsCol name of column with estimated counts
#' @param lengthCol name of column with feature length information
#' @param importer a function used to read in the files
#' @param collatedFiles a character vector of filenames for software which provides
#' abundances and counts in matrix form (e.g. Cufflinks). The files should be, in order,
#' abundances, counts, and a third file with length information
#' @param ignoreTxVersion logical, whether to split the tx id on the '.' character
#' to remove version information, for easier matching with the tx id in gene2tx
#' (default FALSE)
#' @param txi list of matrices of trancript-level abundances, counts, and
#' lengths produced by \code{tximport}, only used by \code{summarizeToGene}
#' 
#' @return a simple list with matrices: abundance, counts, length.
#' A final element 'countsFromAbundance' carries through
#' the character argument used in the tximport call.
#' The length matrix contains the average transcript length for each
#' gene which can be used as an offset for gene-level analysis.
#' Note: tximport does not import bootstrap estimates from kallisto, Salmon, or Sailfish.
#'
#' @describeIn tximport Import tx-level quantifications and summarize
#' abundances, counts and lengths to gene-level (default)
#' or simply output tx-level matrices
#' 
#' @references
#'
#' Charlotte Soneson, Michael I. Love, Mark D. Robinson (2015):
#' Differential analyses for RNA-seq: transcript-level estimates
#' improve gene-level inferences. F1000Research.
#' \url{http://dx.doi.org/10.12688/f1000research.7563.1}
#' 
#' @examples
#'
#' # load data for demonstrating tximport
#' # note that the vignette shows more examples
#' # including how to read in files quickly using the readr package
#' 
#' library(tximportData)
#' dir <- system.file("extdata", package="tximportData")
#' samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
#' files <- file.path(dir,"salmon", samples$run, "quant.sf")
#' names(files) <- paste0("sample",1:6)
#'
#' # tx2gene links transcript IDs to gene IDs for summarization
#' tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))
#'
#' txi <- tximport(files, type="salmon", tx2gene=tx2gene)
#'
#' @importFrom utils read.delim
#' @export
tximport <- function(files,
                     type=c("none","kallisto","salmon","sailfish","rsem"),
                     txIn=TRUE,
                     txOut=FALSE,
                     countsFromAbundance=c("no","scaledTPM","lengthScaledTPM"),
                     tx2gene=NULL,
                     reader=read.delim,
                     geneIdCol,
                     txIdCol,
                     abundanceCol,
                     countsCol,
                     lengthCol,
                     importer,
                     collatedFiles,
                     ignoreTxVersion=FALSE) {

  type <- match.arg(type, c("none","kallisto","salmon","sailfish","rsem"))
  countsFromAbundance <- match.arg(countsFromAbundance, c("no","scaledTPM","lengthScaledTPM"))
  stopifnot(all(file.exists(files)))
  if (!txIn & txOut) stop("txOut only an option when transcript-level data is read in (txIn=TRUE)")
  
  # kallisto presets
  if (type == "kallisto") {
    geneIdCol="gene_id"
    txIdCol <- "target_id"
    abundanceCol <- "tpm"
    countsCol <- "est_counts"
    lengthCol <- "eff_length"
    importer <- reader
    }
  
  # salmon/sailfish presets
  if (type %in% c("salmon","sailfish")) {
    geneIdCol="gene_id"
    txIdCol <- "Name"
    abundanceCol <- "TPM"
    countsCol <- "NumReads"
    lengthCol <- "EffectiveLength"
    importer <- function(x) reader(x, comment='#') 
  }
  
  # rsem presets
  if (type == "rsem") {
    txIn <- FALSE
    geneIdCol <- "gene_id"
    abundanceCol <- "FPKM"
    countsCol <- "expected_count"
    lengthCol <- "effective_length"
    importer <- reader
  }
  
  if (type == "cufflinks") {
    stop("reading from collated files not yet implemented")
  }
  
  # if input is tx-level, need to summarize abundances, counts and lengths to gene-level
  if (txIn) {
    message("reading in files")
    for (i in seq_along(files)) {
      message(i," ",appendLF=FALSE)
      raw <- as.data.frame(importer(files[i]))

      #####################################################################
      # some temporary code for detecting older fishes
      if ((i == 1) &
          (type %in% c("salmon","sailfish")) &
          !("EffectiveLength" %in% names(raw))) {
        lengthCol <- "Length" 
        # because the comment lines have the same comment character
        # as the header, need to name the column names
        importer <- function(x) {
          tmp <- reader(x, comment="#", header=FALSE)
          names(tmp) <- c("Name","Length","TPM","NumReads")
          tmp
        }
        # re-read the first file
        raw <- try(as.data.frame(importer(files[i])), silent=TRUE)
        # if this didn't work, reader is likely read_tsv and
        # different importer() code is needed
        if (inherits(raw, "try-error")) {
          importer <- function(x) {
            reader(x, comment="#", col_names=c("Name","Length","TPM","NumReads"))
          }
          raw <- try(as.data.frame(importer(files[i])))
          if (inherits(raw, "try-error")) stop("tried but couldn't use reader() without error
  user will need to define the importer() as well")
        }
      }
      #####################################################################
      
      # does the table contain gene association or was an external tx2gene table provided?
      if (is.null(tx2gene) & !txOut) {
        # e.g. Cufflinks includes the gene ID in the table
        if (!geneIdCol %in% names(raw)) {
          message()
          stop("

  tximport failed at summarizing to the gene-level.
  Please see 'Solutions' in the Details section of the man page: ?tximport

")
        }
        stopifnot(all(c(lengthCol, abundanceCol) %in% names(raw)))
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
        rownames(mat) <- raw[[txIdCol]]
        colnames(mat) <- names(files)
        abundanceMatTx <- mat
        countsMatTx <- mat
        lengthMatTx <- mat
      }
      abundanceMatTx[,i] <- raw[[abundanceCol]]
      countsMatTx[,i] <- raw[[countsCol]]
      lengthMatTx[,i] <- raw[[lengthCol]]
    }
    message("")

    txi <- list(abundance=abundanceMatTx, counts=countsMatTx, length=lengthMatTx,
                countsFromAbundance="no")

    # if the user requested just the transcript-level data:
    if (txOut) {
      return(txi)
    }

    txi[["countsFromAbundance"]] <- NULL
    txiGene <- summarizeToGene(txi, tx2gene, ignoreTxVersion, countsFromAbundance)
    return(txiGene)  
    
  # e.g. RSEM already has gene-level summaries
  # just combine the gene-level summaries across files
  } else {
  
    message("reading in files")
    for (i in seq_along(files)) {
      message(i," ",appendLF=FALSE)
      raw <- as.data.frame(importer(files[i]))
      stopifnot(all(c(geneIdCol, abundanceCol, lengthCol) %in% names(raw)))
      if (i == 1) {
        mat <- matrix(nrow=nrow(raw),ncol=length(files))
        rownames(mat) <- raw[[geneIdCol]]
        colnames(mat) <- names(files)
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
  return(list(abundance=abundanceMat, counts=countsMat, length=lengthMat,
              countsFromAbundance="no"))
}

# summarizeToGene() splits out the summarization functions
# in tximport(), so it can be called by users to summarize
# transcript-level lists of matrices

#' @describeIn tximport Summarize tx-level matrices to gene-level
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
    if (ntxmissing > 0) message("transcripts missing genes: ", ntxmissing)

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
    countsSum <- colSums(countsMat)
    if (countsFromAbundance == "lengthScaledTPM") {
      newCounts <- abundanceMat * rowMeans(lengthMat)
    } else {
        newCounts <- abundanceMat
      }
    newSum <- colSums(newCounts)
    countsMat <- t(t(newCounts) * (countsSum/newSum))
  }
  
  return(list(abundance=abundanceMat, counts=countsMat, length=lengthMat,
              countsFromAbundance=countsFromAbundance))
}

# function for replacing missing average transcript length values
replaceMissingLength <- function(lengthMat, aveLengthSampGene) {
  nanRows <- which(apply(lengthMat, 1, function(row) any(is.nan(row))))
  if (length(nanRows) > 0) {
    for (i in nanRows) {
      if (all(is.nan(lengthMat[i,]))) {
        # if all samples have 0 abundances for all tx, use the simple average
        lengthMat[i,] <- aveLengthSampGene[i]
      } else {
          # otherwise use the geometric mean of the lengths from the other samples
          idx <- is.nan(lengthMat[i,])
          lengthMat[i,idx] <-  exp(mean(log(lengthMat[i,!idx]), na.rm=TRUE))
        }
    }
  }
  lengthMat
}

# this is much faster than by(), a bit slower than dplyr summarize_each()
## fastby <- function(m, f, fun) {
##   idx <- split(1:nrow(m), f)
##   if (ncol(m) > 1) {
##     t(sapply(idx, function(i) fun(m[i,,drop=FALSE])))
##   } else {
##     matrix(vapply(idx, function(i) fun(m[i,,drop=FALSE], FUN.VALUE=numeric(ncol(m)))),
##            dimnames=list(levels(f), colnames(m)))
##   }
## }

