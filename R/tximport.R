#' Tximport package: import transcript-level quantification data
#' 
#' The tximport package is designed to simplify import of transcript-level
#' abundances (TPM), estimated counts, and effective lengths from
#' a variety of upstream tools, for downstream transcript-level or
#' gene-level analysis. It has no dependencies beyond R, so as to
#' minimize requirements for downstream packages making use of
#' tximport.
#'
#' The main function has the same name as the package:
#' \itemize{
#' \item \code{\link{tximport}} - with key arguments: \code{files}, \code{type}, \code{txOut}, and \code{tx2gene}
#' }
#' 
#' All software-related questions should be posted to the Bioconductor Support Site:
#' 
#' \url{https://support.bioconductor.org}
#'
#' The code can be viewed at the GitHub repository,
#' which also lists the contributor code of conduct:
#'
#' \url{https://github.com/mikelove/tximport}
#' 
#' @references
#'
#' Charlotte Soneson, Michael I. Love, Mark D. Robinson (2015)
#' Differential analyses for RNA-seq: transcript-level estimates
#' improve gene-level inferences. F1000Research.
#' \url{http://doi.org/10.12688/f1000research.7563}
#'
#' @author Charlotte Soneson, Michael I. Love, Mark D. Robinson
#' 
#' @docType package
#' @name tximport-package
#' @aliases tximport-package
#' @keywords package
NULL

#' Import transcript-level abundances and counts
#' for transcript- and gene-level analysis packages
#'
#' \code{tximport} imports transcript-level estimates from various
#' external software and optionally summarizes abundances, counts,
#' and transcript lengths
#' to the gene-level (default) or outputs transcript-level matrices
#' (see \code{txOut} argument).
#'
#' \strong{Inferential replicates:}
#' \code{tximport} will also load in information about inferential replicates --
#' a list of matrices of the Gibbs samples from the posterior, or bootstrap replicates,
#' per sample -- if these data are available in the expected locations relative
#' to the \code{files}.
#' The inferential replicates, stored in \code{infReps} in the output list,
#' are on estimated counts, and therefore follow \code{counts} in the output list.
#' By setting \code{varReduce=TRUE}, the inferential replicate matrices
#' will be replaced by a single matrix with the sample variance per transcript/gene
#' and per sample.
#'
#' \strong{summarizeToGene:}
#' While \code{tximport} summarizes to the gene-level by default, 
#' the user can also perform the import and summarization steps manually,
#' by specifing \code{txOut=TRUE} and then using the function \code{summarizeToGene}.
#' Note however that this is equivalent to \code{tximport} with
#' \code{txOut=FALSE} (the default).
#'
#' \strong{Solutions on summarization:} regarding \code{"tximport failed at summarizing to the gene-level"}:
#'
#' \enumerate{
#'   \item provide a \code{tx2gene} data.frame linking transcripts to genes (more below)
#'   \item avoid gene-level summarization by specifying \code{txOut=TRUE}
#' }
#' 
#' See \code{vignette('tximport')} for example code for generating a
#' \code{tx2gene} data.frame from a \code{TxDb} object.
#' The \code{tx2gene} data.frame should exactly match and be derived from
#' the same set of transcripts used for quantifying (the set of transcript
#' used to create the transcriptome index).
#'
#' \strong{Tximeta:}
#' One automated solution for Salmon or alevin quantification data is to use the
#' \code{tximeta} function in the tximeta Bioconductor package
#' which builds upon and extends \code{tximport}; this solution should
#' work out-of-the-box for human and mouse transcriptomes downloaded
#' from GENCODE, Ensembl, or RefSeq. For other cases, the user
#' should create the \code{tx2gene} manually as shown in the tximport
#' vignette.
#'
#' \strong{On tx2gene construction:}
#' Note that the \code{keys} and \code{select} functions used
#' to create the \code{tx2gene} object are documented
#' in the man page for \link[AnnotationDbi]{AnnotationDb-class} objects
#' in the AnnotationDbi package (TxDb inherits from AnnotationDb).
#' For further details on generating TxDb objects from various inputs
#' see \code{vignette('GenomicFeatures')} from the GenomicFeatures package.
#'
#' \strong{alevin:}
#' The \code{alevinArgs} argument includes some alevin-specific arguments.
#' This optional argument is a list with any or all of the following named logical variables:
#' \code{filterBarcodes}, \code{tierImport}, and \code{forceSlow}.
#' The variables are described as follows (with default values in parens):
#' \code{filterBarcodes} (FALSE) import only cell barcodes listed in
#' \code{whitelist.txt};
#' \code{tierImport} (FALSE) import the tier information in addition to counts;
#' \code{forceSlow} (FALSE) force the use of the slower import R code
#' even if \code{eds} is installed;
#' \code{dropMeanVar} (FALSE) don't import inferential mean and variance
#' matrices even if they exist (also skips inferential replicates)
#' For \code{type="alevin"} all arguments other than \code{files},
#' \code{dropInfReps}, and \code{alevinArgs} are ignored.
#' Note that \code{files} should point to a single \code{quants_mat.gz} file,
#' in the directory structure created by the alevin software
#' (e.g. do not move the file or delete the other important files).
#' Note that importing alevin quantifications will be much faster by first
#' installing the \code{eds} package, which contains a C++ importer
#' for alevin's EDS format.
#' For alevin, \code{tximport} is importing the gene-by-cell matrix of counts,
#' as \code{txi$counts}, and effective lengths are not estimated.
#' \code{txi$mean} and \code{txi$variance} may also be imported if
#' inferential replicates were used, as well as inferential replicates
#' if these were output by alevin.
#' Length correction should not be applied to datasets where there
#' is not an expected correlation of counts and feature length.
#' 
#' @param files a character vector of filenames for the transcript-level abundances
#' @param type character, the type of software used to generate the abundances.
#' Options are "salmon", "sailfish", "alevin", "kallisto", "rsem", "stringtie", or "none".
#' This argument is used to autofill the arguments below (geneIdCol, etc.)
#' "none" means that the user will specify these columns. Be aware that
#' specifying \code{type} other than "none" will ignore the arguments below
#' (geneIdCol, etc.)
#' @param txIn logical, whether the incoming files are transcript level (default TRUE)
#' @param txOut logical, whether the function should just output
#' transcript-level (default FALSE)
#' @param countsFromAbundance character, either "no" (default), "scaledTPM",
#' "lengthScaledTPM", or "dtuScaledTPM".
#' Whether to generate estimated counts using abundance estimates:
#' \itemize{
#'   \item scaled up to library size (scaledTPM),
#'   \item scaled using the average transcript length over samples
#'         and then the library size (lengthScaledTPM), or
#'   \item scaled using the median transcript length among isoforms of a gene,
#'         and then the library size (dtuScaledTPM). 
#' }
#' dtuScaledTPM is designed for DTU analysis in combination with \code{txOut=TRUE},
#' and it requires specifing a \code{tx2gene} data.frame.
#' dtuScaledTPM works such that within a gene, values from all samples and
#' all transcripts get scaled by the same fixed median transcript length.
#' If using scaledTPM, lengthScaledTPM, or geneLengthScaledTPM, 
#' the counts are no longer correlated across samples with transcript length,
#' and so the length offset matrix should not be used.
#' @param tx2gene a two-column data.frame linking transcript id (column 1)
#' to gene id (column 2).
#' the column names are not relevant, but this column order must be used. 
#' this argument is required for gene-level summarization, and the tximport
#' vignette describes how to construct this data.frame (see Details below).
#' An automated solution to avoid having to create \code{tx2gene} if
#' one has quantified with Salmon or alevin with human or mouse transcriptomes
#' is to use the \code{tximeta} function from the tximeta Bioconductor package.
#' @param varReduce whether to reduce per-sample inferential replicates
#' information into a matrix of sample variances \code{variance} (default FALSE).
#' alevin computes inferential variance by default for bootstrap
#' inferential replicates, so this argument is ignored/not necessary
#' @param dropInfReps whether to skip reading in inferential replicates
#' (default FALSE). For alevin, \code{tximport} will still read in the
#' inferential variance matrix if it exists
#' @param infRepStat a function to re-compute counts and abundances from the
#' inferential replicates, e.g. \code{matrixStats::rowMedians} to re-compute counts 
#' as the median of the inferential replicates. The order of operations is:
#' first counts are re-computed, then abundances are re-computed.
#' Following this, if \code{countsFromAbundance} is not "no",
#' \code{tximport} will again re-compute counts from the re-computed abundances.
#' \code{infRepStat} should operate on rows of a matrix. (default is NULL)
#' @param ignoreTxVersion logical, whether to split the tx id on the '.' character
#' to remove version information to facilitate matching with the tx id in \code{tx2gene}
#' (default FALSE)
#' @param ignoreAfterBar logical, whether to split the tx id on the '|' character
#' to facilitate matching with the tx id in \code{tx2gene} (default FALSE).
#' if \code{txOut=TRUE} it will strip the text after '|' on the rownames
#' of the matrices
#' @param geneIdCol name of column with gene id. if missing, the \code{tx2gene}
#' argument can be used. Note that this argument and the other four "...Col"
#' arguments below are ignored unless \code{type="none"} 
#' @param txIdCol name of column with tx id
#' @param abundanceCol name of column with abundances (e.g. TPM or FPKM)
#' @param countsCol name of column with estimated counts
#' @param lengthCol name of column with feature length information
#' @param importer a function used to read in the files
#' @param existenceOptional logical, should tximport not check if files exist before attempting
#' import (default FALSE, meaning files must exist according to \code{file.exists})
#' @param sparse logical, whether to try to import data sparsely (default is FALSE).
#' Initial implementation for \code{txOut=TRUE}, \code{countsFromAbundance="no"}
#' or \code{"scaledTPM"}, no inferential replicates. Only counts matrix
#' is returned (and abundance matrix if using \code{"scaledTPM"})
#' @param sparseThreshold the minimum threshold for including a count as a
#' non-zero count during sparse import (default is 1)
#' @param readLength numeric, the read length used to calculate counts from
#' StringTie's output of coverage. Default value (from StringTie) is 75.
#' The formula used to calculate counts is:
#' \code{cov * transcript length / read length}
#' @param alevinArgs named list, with logical elements \code{filterBarcodes},
#' \code{tierImport}, \code{forceSlow}, \code{dropMeanVar}.
#' See Details for definitions.
#' 
#' @return A simple list containing matrices: abundance, counts, length.
#' Another list element 'countsFromAbundance' carries through
#' the character argument used in the tximport call.
#' The length matrix contains the average transcript length for each
#' gene which can be used as an offset for gene-level analysis.
#' If detected, and \code{txOut=TRUE}, inferential replicates for
#' each sample will be imported and stored as a list of matrices,
#' itself an element \code{infReps} in the returned list.
#' An exception is alevin, in which the \code{infReps} are a list
#' of bootstrap replicate matrices, where each matrix has
#' genes as rows and cells as columns.
#' If \code{varReduce=TRUE} the inferential replicates will be summarized
#' according to the sample variance, and stored as a matrix \code{variance}.
#' alevin already computes the variance of the bootstrap inferential replicates
#' and so this is imported without needing to specify \code{varReduce=TRUE}.
#' 
#' @references
#'
#' Charlotte Soneson, Michael I. Love, Mark D. Robinson (2015)
#' Differential analyses for RNA-seq: transcript-level estimates
#' improve gene-level inferences. F1000Research.
#' \url{http://doi.org/10.12688/f1000research.7563}
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
#' files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
#' names(files) <- paste0("sample",1:6)
#'
#' # tx2gene links transcript IDs to gene IDs for summarization
#' tx2gene <- read.csv(file.path(dir, "tx2gene.gencode.v27.csv"))
#'
#' txi <- tximport(files, type="salmon", tx2gene=tx2gene)
#'
#' @importFrom utils read.delim capture.output head compareVersion packageVersion
#' @importFrom stats median
#' @importFrom methods is
#'
#' @export
tximport <- function(files,
                     type=c("none","salmon","sailfish","alevin","kallisto","rsem","stringtie", "piscem"),
                     txIn=TRUE,
                     txOut=FALSE,
                     countsFromAbundance=c("no","scaledTPM","lengthScaledTPM","dtuScaledTPM"),
                     tx2gene=NULL,
                     varReduce=FALSE,
                     dropInfReps=FALSE,
                     infRepStat=NULL,
                     ignoreTxVersion=FALSE,
                     ignoreAfterBar=FALSE,
                     geneIdCol,
                     txIdCol,
                     abundanceCol,
                     countsCol,
                     lengthCol,
                     importer=NULL,
                     existenceOptional=FALSE,
                     sparse=FALSE,
                     sparseThreshold=1,
                     readLength=75,
                     alevinArgs=NULL) {
  
  # inferential replicate importer
  infRepImporter <- NULL
  
  type <- match.arg(type)
  countsFromAbundance <- match.arg(countsFromAbundance)
  if (countsFromAbundance == "dtuScaledTPM") {
    stopifnot(txOut)
    if (is.null(tx2gene)) stop("'dtuScaledTPM' requires 'tx2gene' input")
  }
  
  if (!existenceOptional) stopifnot(all(file.exists(files)))
  if (!txIn & txOut) stop("txOut only an option when transcript-level data is read in (txIn=TRUE)")
  
  stopifnot(length(files) > 0)
  kallisto.h5 <- basename(files[1]) == "abundance.h5"
  if (type == "kallisto" & !kallisto.h5) {
    message("Note: importing `abundance.h5` is typically faster than `abundance.tsv`")
  }
  
  if (type=="rsem" & txIn & grepl("genes", files[1])) {
    message("It looks like you are importing RSEM genes.results files, setting txIn=FALSE")
    txIn <- FALSE
  }
  
  # in order to use infRepStat, we need inf reps, and doesn't work with alevin or sparse code
  if (!is.null(infRepStat)) {
    if (dropInfReps) stop("infRepStat requires infReps")
    if (type == "alevin") stop("infRepStat does not currently work with alevin input")
    if (sparse) stop("infRepStat does not currently work with sparse output")
  }
  
  # special alevin code
  if (type=="alevin") {
    # unpacking alevinArgs
    if (is.null(alevinArgs)) {
      alevinArgs <- list(filterBarcodes=FALSE, tierImport=FALSE,
                         forceSlow=FALSE, dropMeanVar=FALSE)
    }
    stopifnot(is(alevinArgs, "list"))
    stopifnot(all(sapply(alevinArgs, is.logical)))
    alevinArgNms <- c("filterBarcodes","tierImport","forceSlow","dropMeanVar")
    stopifnot(all(names(alevinArgs) %in% alevinArgNms))
    filterBarcodes <- if (is.null(alevinArgs$filterBarcodes)) FALSE else alevinArgs$filterBarcodes
    tierImport <- if (is.null(alevinArgs$tierImport)) FALSE else alevinArgs$tierImport
    forceSlow <- if (is.null(alevinArgs$forceSlow)) FALSE else alevinArgs$forceSlow
    dropMeanVar <- if (is.null(alevinArgs$dropMeanVar)) FALSE else alevinArgs$dropMeanVar    
    if (length(files) > 1) stop("alevin import currently only supports a single experiment")
    # check that alevin is >= version 0.14.0
    if (compareVersion(getAlevinVersion(files), "0.14.0") == -1) {
      stop("use of tximport version >= 1.18 requires alevin version >= 0.14")
    }
    mat <- readAlevin(files, dropInfReps, filterBarcodes, tierImport, forceSlow, dropMeanVar)
    # if 'mat' is not a list, it is a matrix with the counts
    if (!is.list(mat)) {
      txi <- list(abundance=NULL, counts=mat)
    } else {
      # otherwise, 'mat' is a list with various information to pass on
      txi <- list(abundance=NULL, counts=mat$counts, mean=mat$mean, variance=mat$variance)
      if (tierImport) {
        txi$tier <- mat$tier
      }
      if ("infReps" %in% names(mat)) {
        txi$infReps <- mat$infReps
      }
    }
    txi$length <- NULL
    txi$countsFromAbundance="no"
    return(txi)
  }
  
  readrStatus <- FALSE
  if (is.null(importer) & !kallisto.h5) {
    if (!requireNamespace("readr", quietly=TRUE)) {
      message("reading in files with read.delim (install 'readr' package for speed up)")
      importer <- read.delim
    } else {
      # this fix contributed by @ATpoint
      # readr won't work on very stripped down machines, e.g. Docker images w/o timezone
      timeZoneMissing <- length(suppressWarnings(OlsonNames())) == 0 | is.na(Sys.timezone())
      if (timeZoneMissing) {
        message("reading in files with read.delim ('readr' installed but won't work w/o timezones)")
        importer <- read.delim
      } else {
        message("reading in files with read_tsv")
        readrStatus <- TRUE
      }
    }
  }
  
  # salmon/sailfish presets
  if (type %in% c("salmon","sailfish")) {
    txIdCol <- "Name"
    abundanceCol <- "TPM"
    countsCol <- "NumReads"
    lengthCol <- "EffectiveLength"
    if (readrStatus & is.null(importer)) {
      col.types <- readr::cols(
        readr::col_character(),readr::col_integer(),readr::col_double(),readr::col_double(),readr::col_double()
      )
      importer <- function(x) readr::read_tsv(x, progress=FALSE, col_types=col.types)
    }
    infRepImporter <- if (dropInfReps) { NULL } else { function(x) readInfRepFish(x, type) }
  }

  # piscem presets
  if (type == "piscem") {
    txIdCol <- "target_name"
    lengthCol <- "eelen"
    abundanceCol <- "tpm"
    countsCol <- "ecount"
    col.types <- readr::cols(
                             readr::col_character(),readr::col_integer(),readr::col_double(),readr::col_double(),readr::col_double()
    )
    importer <- function(x) readr::read_tsv(x, progress=FALSE, col_types=col.types)
    infRepImporter <- if (dropInfReps) { NULL } else { readInfRepPiscem }
  }
  
  # kallisto presets
  if (type == "kallisto") {
    txIdCol <- "target_id"
    abundanceCol <- "tpm"
    countsCol <- "est_counts"
    lengthCol <- "eff_length"
    if (kallisto.h5) {
      importer <- read_kallisto_h5
    } else if (readrStatus & is.null(importer)) {
      col.types <- readr::cols(
        readr::col_character(),readr::col_integer(),readr::col_double(),readr::col_double(),readr::col_double()
      )
      importer <- function(x) readr::read_tsv(x, progress=FALSE, col_types=col.types)
    }
    infRepImporter <- if (dropInfReps) { NULL } else { readInfRepKallisto }
  }
  
  # rsem presets
  if (type == "rsem") {
    if (txIn) {
      txIdCol <- "transcript_id"
      abundanceCol <- "TPM"
      countsCol <- "expected_count"
      lengthCol <- "effective_length"
      if (readrStatus & is.null(importer)) {
        col.types <- readr::cols(
          readr::col_character(),readr::col_character(),readr::col_integer(),readr::col_double(),
          readr::col_double(),readr::col_double(),readr::col_double(),readr::col_double()
        )
        importer <- function(x) readr::read_tsv(x, progress=FALSE, col_types=col.types)
      }
    } else {
      geneIdCol <- "gene_id"
      abundanceCol <- "TPM"
      countsCol <- "expected_count"
      lengthCol <- "effective_length"
      if (readrStatus & is.null(importer)) {
        col.types <- readr::cols(
          readr::col_character(),readr::col_character(),readr::col_double(),readr::col_double(),
          readr::col_double(),readr::col_double(),readr::col_double()
        )
        importer <- function(x) readr::read_tsv(x, progress=FALSE, col_types=col.types)
      }
    }
  }
  
  if (type == c("stringtie")) {
    txIdCol <- "t_name"
    geneIdCol <- "gene_name"
    abundanceCol <- "FPKM"
    countsCol <- "cov"
    lengthCol <- "length"
    if (readrStatus & is.null(importer)) {
      col.types <- readr::cols(
        readr::col_character(),readr::col_character(),readr::col_character(),readr::col_integer(),readr::col_integer(),readr::col_character(),readr::col_integer(),readr::col_integer(),readr::col_character(),readr::col_character(),readr::col_double(),readr::col_double()
      )
      importer <- function(x) readr::read_tsv(x, progress=FALSE, col_types=col.types)
    }
  }
  
  infRepType <- "none"
  if (type %in% c("salmon", "sailfish", "piscem", "kallisto") & !dropInfReps) {
    # if summarizing to gene-level, need the full matrices passed to summarizeToGene
    infRepType <- if (varReduce & txOut) { "var" } else { "full" }
  }
  
  # special code for RSEM gene.results files.
  # RSEM gene-level is the only case of !txIn
  if (!txIn) {
    txi <- computeRsemGeneLevel(files, importer, geneIdCol, abundanceCol, countsCol, lengthCol, countsFromAbundance)
    message("")
    return(txi)
  }
  
  # if external tx2gene table not provided, send user to vignette
  if (is.null(tx2gene) & !txOut) {
    summarizeFail() # ...long message in helper.R
  }
  
  # trial run of inferential replicate info
  repInfo <- NULL
  if (infRepType != "none") {
    repInfo <- if (type == "piscem") { infRepImporter(files[1]) } else { infRepImporter(dirname(files[1])) }
    # if we didn't find inferential replicate info
    if (is.null(repInfo)) {
      infRepType <- "none"
    }
  }
  
  if (sparse) {
    if (!requireNamespace("Matrix", quietly=TRUE)) {
      stop("sparse import requires core R package `Matrix`")
    }
    message("importing sparsely, only counts and abundances returned, support limited to
txOut=TRUE, CFA either 'no' or 'scaledTPM', and no inferential replicates")
    stopifnot(txOut)
    stopifnot(infRepType == "none")
    stopifnot(countsFromAbundance %in% c("no","scaledTPM"))
  }
  
  ######################################################
  # the rest of the code assumes transcript-level input:
  
  ### --- BEGIN --- loop over files reading in columns / inf reps ###
  for (i in seq_along(files)) {
    message(i," ",appendLF=FALSE)
    # import and convert quantification info to data.frame
    raw <- as.data.frame(importer(files[i]))
    # import inferential replicate info
    if (infRepType != "none") {
      repInfo <- if (type == "piscem") { infRepImporter(files[i]) } else { infRepImporter(dirname(files[i])) }
    } else {
      repInfo <- NULL
    }
    # check for columns
    stopifnot(all(c(abundanceCol, countsCol, lengthCol) %in% names(raw)))
    # check for same-across-samples
    if (i == 1) {
      txId <- raw[[txIdCol]]
    } else {
      stopifnot(all(txId == raw[[txIdCol]]))
    }
    # if importing dense matrices
    if (!sparse) {
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
      # if infRepStat was specified, re-compute counts and abundances
      if (!is.null(infRepStat)) {
        countsMatTx[,i] <- infRepStat(repInfo$reps)
        tpm <- countsMatTx[,i] / lengthMatTx[,i]
        abundanceMatTx[,i] <- tpm * 1e6 / sum(tpm)
      }
    } else {
      # try importing sparsely
      if (i == 1) {     
        txId <- raw[[txIdCol]]
        countsListI <- list()
        countsListX <- list()
        abundanceListX <- list()
        numNonzero <- c()
      }
      stopifnot(all(txId == raw[[txIdCol]]))
      sparse.idx <- which(raw[[countsCol]] >= sparseThreshold)
      countsListI <- c(countsListI, sparse.idx)
      countsListX <- c(countsListX, raw[[countsCol]][sparse.idx])
      numNonzero <- c(numNonzero, length(sparse.idx))
      if (countsFromAbundance == "scaledTPM") {
        abundanceListX <- c(abundanceListX, raw[[abundanceCol]][sparse.idx])
      }
    }
  }
  ### --- END --- loop over files ###
  
  # compile sparse matrices
  if (sparse) {
    countsMatTx <- Matrix::sparseMatrix(i=unlist(countsListI),
                                        j=rep(seq_along(numNonzero), numNonzero),
                                        x=unlist(countsListX),
                                        dims=c(length(txId), length(files)),
                                        dimnames=list(txId, names(files)))
    if (countsFromAbundance == "scaledTPM") {
      abundanceMatTx <- Matrix::sparseMatrix(i=unlist(countsListI),
                                             j=rep(seq_along(numNonzero), numNonzero),
                                             x=unlist(abundanceListX),
                                             dims=c(length(txId), length(files)),
                                             dimnames=list(txId, names(files)))
    } else {
      abundanceMatTx <- NULL
    }
    lengthMatTx <- NULL
  }
  
  # propagate names to inferential replicate list
  if (infRepType == "full") {
    if (length(infRepMatTx) != length(files)) {
      stop("Note: not all samples contain inferential replicates.
  tximport can only import data when either all or no samples
  contain inferential replicates. Instead first subset to the
  set of samples that all contain inferential replicates.")
    }
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
  
  # stringtie outputs coverage, here we turn into counts
  if (type == "stringtie") {
    # here "counts" is still just coverage, this formula gives back original counts
    txi$counts <- txi$counts * txi$length / readLength
  }
  
  if (type == "rsem") {
    # protect against 0 bp length transcripts
    txi$length[txi$length < 1] <- 1
  }
  
  # two main outputs, based on choice of txOut:
  
  # 1) if the user requested just the transcript-level data, return it now
  if (txOut) {
    # if countsFromAbundance in {scaledTPM, lengthScaledTPM, or dtuScaledTPM}
    if (countsFromAbundance != "no") {
      # for dtuScaledTPM, pretend we're doing lengthScaledTPM w/ an altered length matrix.
      # note that we will still output txi$countsFromAbundance set to "dtuScaledTPM"
      length4CFA <- txi$length # intermediate version of the length matrix
      if (countsFromAbundance == "dtuScaledTPM") {
        length4CFA <- medianLengthOverIsoform(length4CFA, tx2gene, ignoreTxVersion, ignoreAfterBar)
        countsFromAbundance <- "lengthScaledTPM" 
      }
      # function for computing all 3 countsFromAbundance methods:
      txi$counts <- makeCountsFromAbundance(countsMat=txi$counts,
                                            abundanceMat=txi$abundance,
                                            lengthMat=length4CFA,
                                            countsFromAbundance=countsFromAbundance)
    }
    # request from Issue 40: remove text after '|' for txOut=TRUE
    if (ignoreAfterBar) {
      for (matNm in c("counts","abundance","length")) {
        rowNms <- rownames(txi[[matNm]])
        rownames(txi[[matNm]]) <- sub("\\|.*", "", rowNms)
      }
    }
    return(txi)
  }
  
  # 2) otherwise, summarize to the gene-level
  txi[["countsFromAbundance"]] <- NULL
  txiGene <- summarizeToGene(txi, tx2gene, varReduce, ignoreTxVersion, ignoreAfterBar, countsFromAbundance)
  return(txiGene)
  
}


# split out this special code for RSEM with gene-level input
# (all other input is transcript-level)
computeRsemGeneLevel <- function(files, importer,
                                 geneIdCol, abundanceCol, countsCol, lengthCol,
                                 countsFromAbundance) {
  # RSEM already has gene-level summaries
  # so we just combine the gene-level summaries across files
  if (countsFromAbundance != "no") {
    warning("countsFromAbundance other than 'no' requires transcript-level estimates")
  }
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
  txi <- list(abundance=abundanceMat, counts=countsMat, length=lengthMat,
              countsFromAbundance="no")
  return(txi)
}
