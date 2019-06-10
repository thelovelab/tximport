# from http://stackoverflow.com/questions/25099825/row-wise-variance-of-a-matrix-in-r
rowVars <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

# Read inferential replicate information from salmon / sailfish
#
# SF ver >= 0.9.0, Salmon ver >= 0.8.0
#
# fish_dir = path to a sailfish output directory
readInfRepFish <- function(fish_dir, meth) {
  # aux_info is the default auxiliary directory in salmon
  # aux is the default directory in sailfish
  aux_dir <- "aux_info"

  if (!requireNamespace("jsonlite", quietly=TRUE)) {
    stop("importing inferential replicates for Salmon or Sailfish requires package `jsonlite`.
  to skip this step, set dropInfReps=TRUE")
  }
  
  # if the default is overwritten, then use that instead
  jsonPath <- file.path(fish_dir, "cmd_info.json")
  if (!file.exists(jsonPath)) {
    return(NULL)
  }
  
  cmd_info <- jsonlite::fromJSON(jsonPath)
  if (is.element("auxDir", names(cmd_info))) {
    aux_dir <- cmd_info$auxDir
  }
  auxPath <- file.path(fish_dir, aux_dir)
  if (!file.exists(auxPath)) {
    return(NULL)
  }

  # get all of the meta info
  minfo <- jsonlite::fromJSON(file.path(auxPath, "meta_info.json"))

  if ("salmon_version" %in% names(minfo)) {
    stopifnot(package_version(minfo$salmon_version) >= "0.8.0")
  }
  if ("sailfish_version" %in% names(minfo)) {
    stopifnot(package_version(minfo$sailfish_version) >= "0.9.0")
  }
  
  sampType <- NULL
  # check if we have explicitly recorded the type of posterior sample
  # (salmon >= 0.7.3)
  if (is.element("samp_type", names(minfo))) {
    sampType <- minfo$samp_type
  }

  # load bootstrap data if it exists
  knownSampleTypes <- c("gibbs", "bootstrap")
  numBoot <- minfo$num_bootstraps
  if (numBoot > 0) {
    bootCon <- gzcon(file(file.path(auxPath, 'bootstrap', 'bootstraps.gz'), "rb"))
    ##
    # Gibbs samples *used* to be integers, and bootstraps were doubles
    # Now, however, both types of samples are doubles.  The code below 
    # tries to load doubles first, but falls back to integers if it fails.
    ##
    if("num_valid_targets" %in% names(minfo)) {
      minfo$num_targets = minfo$num_valid_targets
    }
    expected.n <- minfo$num_targets * minfo$num_bootstraps
    boots <- tryCatch({
      bootsIn <- readBin(bootCon, "double", n = expected.n)
      stopifnot(length(bootsIn) == expected.n)
      bootsIn
    }, error=function(...) {
      # close and re-open the connection to reset the file
      close(bootCon)
      bootCon <- gzcon(file(file.path(auxPath, 'bootstrap', 'bootstraps.gz'), "rb"))
      readBin(bootCon, "integer", n = expected.n)
    })
    close(bootCon)

    # rows are transcripts, columns are bootstraps
    dim(boots) <- c(minfo$num_targets, minfo$num_bootstraps)
    vars <- rowVars(boots)
    return(list(vars=vars, reps=boots))
  } else {
    return(NULL)
  }
}

readInfRepKallisto <- function(bear_dir) {
  h5File <- file.path(bear_dir, "abundance.h5")
  if (!file.exists(h5File)) return(NULL)
  if (!requireNamespace("rhdf5", quietly=TRUE)) {
    stop("reading kallisto results from hdf5 files requires Bioconductor package `rhdf5`")
  }
  groups <- rhdf5::h5ls(h5File)
  numBoot <- length(groups$group[groups$group == "/bootstrap"])
  if (numBoot > 0) {
    boots <- rhdf5::h5read(file.path(bear_dir, "abundance.h5"), "/bootstrap")
    numTxp <- length(boots$bs0)
    bootMat <- matrix(nrow=numTxp, ncol=numBoot)
    for (bsn in seq_len(numBoot)) {
      bootMat[,bsn] <- boots[bsn][[1]]
    }
    vars <- rowVars(bootMat)
    return(list(vars=vars, reps=bootMat))
  } else {
    return(NULL)
  }
}
