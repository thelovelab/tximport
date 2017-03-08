#' Import transcript-level abundances and counts
#' for transcript- and gene-level analysis packages
#'
#' \code{tximport} imports transcript-level estimates from various
#' external software and optionally summarizes abundances, counts,
#' and transcript lengths
#' to the gene-level (default) or outputs transcript-level matrices
#' (see \code{txOut} argument).
#'
#' \code{tximport} will also load in information about inferential replicates --
#' a list of matrices of the Gibbs samples from the posterior, or bootstrap replicates,
#' per sample -- if these data are available in the expected locations
#' relative to the \code{files}, and if \code{txOut=TRUE}.
#' The inferential replicates, stored in \code{infReps} in the output list,
#' are on estimated counts, and therefore follow \code{counts} in the output list.
#' By setting \code{varReduce=TRUE}, the inferential replicate matrices
#' will be replaced by a single matrix with the sample variance
#' per transcript and per sample.
#' Inferential replicate information is not summarized to the gene level.
#'
#' While \code{tximport} summarizes to the gene-level by default, 
#' the user can also perform the import and summarization steps manually,
#' by specifing \code{txOut=TRUE} and then using the function
#' \code{summarizeToGene}.
#' Note however that this is equivalent to \code{tximport} with
#' \code{txOut=FALSE} (the default).
#'
#' Solutions to the error "tximport failed at summarizing to the gene-level":
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
#' @param files a character vector of filenames for the transcript-level abundances
#' @param type character, the type of software used to generate the abundances.
#' Options are "salmon", "sailfish", "kallisto", "rsem".
#' This argument is used to autofill the arguments below (geneIdCol, etc.)
#' "none" means that the user will specify these columns.
#' @param txIn logical, whether the incoming files are transcript level (default TRUE)
#' @param txOut logical, whether the function should just output
#' transcript-level (default FALSE)
#' @param countsFromAbundance character, either "no" (default), "scaledTPM",
#' or "lengthScaledTPM",
#' for whether to generate estimated counts using abundance estimates
#' scaled up to library size (scaledTPM) or additionally scaled using
#' the average transcript length over samples and
#' the library size (lengthScaledTPM). if using scaledTPM or lengthScaledTPM, 
#' then the counts are no longer correlated with average transcript length,
#' and so the length offset matrix should not be used.
#' @param tx2gene a two-column data.frame linking transcript id (column 1)
#' to gene id (column 2).
#' the column names are not relevant, but this column order must be used. 
#' this argument is required for gene-level summarization for methods
#' that provides transcript-level estimates only
#' (kallisto, Salmon, Sailfish)
#' @param varReduce whether to reduce per-sample inferential replicates
#' information into a matrix of sample variances \code{variance} (default FALSE)
#' @param dropInfReps whether to skip reading in inferential replicates
#' (default FALSE)
#' @param ignoreTxVersion logical, whether to split the tx id on the '.' character
#' to remove version information, for easier matching with the tx id in gene2tx
#' (default FALSE)
#' @param geneIdCol name of column with gene id. if missing,
#' the gene2tx argument can be used
#' @param txIdCol name of column with tx id
#' @param abundanceCol name of column with abundances (e.g. TPM or FPKM)
#' @param countsCol name of column with estimated counts
#' @param lengthCol name of column with feature length information
#' @param importer a function used to read in the files
#' @param txi list of matrices of trancript-level abundances, counts, and
#' lengths produced by \code{tximport}, only used by \code{summarizeToGene}
#' 
#' @return a simple list containing matrices: abundance, counts, length.
#' Another list element 'countsFromAbundance' carries through
#' the character argument used in the tximport call.
#' If detected, and \code{txOut=TRUE}, inferential replicates for
#' each sample will be imported and stored as a list of matrices,
#' itself an element \code{infReps} in the returned list.
#' If \code{varReduce=TRUE} the inferential replicates will be summarized
#' according to the sample variance, and stored as a matrix \code{variance}.
#' The length matrix contains the average transcript length for each
#' gene which can be used as an offset for gene-level analysis.
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
#' @importFrom utils read.delim capture.output
#'
#' @describeIn tximport Import estimates of abundances and counts
#' @export
tximport <- function(files,
                     type=c("none","salmon","sailfish","kallisto","rsem"),
                     txIn=TRUE,
                     txOut=FALSE,
                     countsFromAbundance=c("no","scaledTPM","lengthScaledTPM"),
                     tx2gene=NULL,
                     varReduce=FALSE,
                     dropInfReps=FALSE,
                     ignoreTxVersion=FALSE,
                     geneIdCol,
                     txIdCol,
                     abundanceCol,
                     countsCol,
                     lengthCol,
                     importer=NULL) {

  # inferential replicate importer
  infRepImporter <- NULL

  type <- match.arg(type, c("none","salmon","sailfish","kallisto","rsem"))
  countsFromAbundance <- match.arg(countsFromAbundance, c("no","scaledTPM","lengthScaledTPM"))

  stopifnot(all(file.exists(files)))
  if (!txIn & txOut) stop("txOut only an option when transcript-level data is read in (txIn=TRUE)")

  kallisto.h5 <- basename(files[1]) == "abundance.h5"
  if (type == "kallisto" & !kallisto.h5) {
    message("Note: importing `abundance.h5` is typically faster than `abundance.tsv`")
  }
  
  readrStatus <- FALSE
  if (is.null(importer) & !kallisto.h5) {
    if (!requireNamespace("readr", quietly=TRUE)) {
      message("reading in files with read.delim (install 'readr' package for speed up)")
      importer <- read.delim
    } else {
      message("reading in files with read_tsv")
      readrStatus <- TRUE
      importer <- function(x) readr::read_tsv(x, progress=FALSE, col_types=readr::cols())
    }
  }
  
  # salmon/sailfish presets
  if (type %in% c("salmon","sailfish")) {
    geneIdCol="gene_id"
    txIdCol <- "Name"
    abundanceCol <- "TPM"
    countsCol <- "NumReads"
    lengthCol <- "EffectiveLength"
    infRepImporter <- if (dropInfReps) { NULL } else { function(x) readInfRepFish(x, type) }
  }

  # kallisto presets
  if (type == "kallisto") {
    geneIdCol="gene_id"
    txIdCol <- "target_id"
    abundanceCol <- "tpm"
    countsCol <- "est_counts"
    lengthCol <- "eff_length"
    if (kallisto.h5) {
      importer <- read_kallisto_h5
    }
    infRepImporter <- if (dropInfReps) { NULL } else { readInfRepKallisto }
  }
  
  # rsem presets
  if (type == "rsem") {
    txIn <- FALSE
    geneIdCol <- "gene_id"
    abundanceCol <- "FPKM"
    countsCol <- "expected_count"
    lengthCol <- "effective_length"
  }

  infRepType <- "none"
  if (type %in% c("salmon", "sailfish", "kallisto") & !dropInfReps) {
    infRepType <- if (varReduce) { "var" } else { "full" }
  }
  
  # if input is tx-level, need to summarize abundances, counts and lengths to gene-level
  if (txIn) {
    for (i in seq_along(files)) {
      message(i," ",appendLF=FALSE)

      raw <- as.data.frame(importer(files[i]))

      # if we expect inferential replicate info
      repInfo <- NULL
      if (infRepType != "none") {
        repInfo <- infRepImporter(dirname(files[i]))
        # if we didn't find inferential replicate info
        if (is.null(repInfo)) {
          infRepType <- "none"
        }
      }

      # if external tx2gene table not provided, send user to vignette
      if (is.null(tx2gene) & !txOut) {
        summarizeFail() # ...long message in helper.R
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
        if (infRepType == "var") {
          varMatTx <- mat
        } else if (infRepType == "full") {
          infRepMatTx <- list()
        }
      }
      abundanceMatTx[,i] <- raw[[abundanceCol]]
      countsMatTx[,i] <- raw[[countsCol]]
      lengthMatTx[,i] <- raw[[lengthCol]]
      if (infRepType == "var") {
        varMatTx[,i] <- repInfo$vars
      } else if (infRepType == "full") {
        infRepMatTx[[i]] <- repInfo$reps
      }
    }

    # propagate names to inferential replicate list
    if (infRepType == "full") {
      names(infRepMatTx) <- names(files)
    }
    
    message("")

    # if there is no information about inferential replicates
    if (infRepType == "none") {
      txi <- list(abundance=abundanceMatTx,
                  counts=countsMatTx,
                  length=lengthMatTx,
                  countsFromAbundance=countsFromAbundance)
    } else if (infRepType == "var") {
    # if we're keeping only the variance from inferential replicates
      txi <- list(abundance=abundanceMatTx,
                  counts=countsMatTx,
                  variance=varMatTx,
                  length=lengthMatTx,
                  countsFromAbundance=countsFromAbundance)
    } else if (infRepType == "full") {
    # if we're keeping the full samples from inferential replicates
      txi <- list(abundance=abundanceMatTx,
                  counts=countsMatTx,
                  infReps=infRepMatTx,
                  length=lengthMatTx,
                  countsFromAbundance=countsFromAbundance)
    }

    # if the user requested just the transcript-level data:
    if (txOut) {
      if (countsFromAbundance != "no") {
        txi$counts <- makeCountsFromAbundance(txi$counts,
                                              txi$abundance,
                                              txi$length,
                                              countsFromAbundance)
      }
      return(txi)
    }

    txi[["countsFromAbundance"]] <- NULL
    txiGene <- summarizeToGene(txi, tx2gene, ignoreTxVersion, countsFromAbundance)
    return(txiGene)  
    
  # e.g. RSEM already has gene-level summaries
  # just combine the gene-level summaries across files
  } else {
  
    for (i in seq_along(files)) {
      message(i," ",appendLF=FALSE)

      out <- capture.output({
        raw <- as.data.frame(importer(files[i]))
      }, type="message")
      
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

