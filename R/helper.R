# function for generating counts from abundances
makeCountsFromAbundance <- function(countsMat, abundanceMat, lengthMat, countsFromAbundance) {
  countsSum <- colSums(countsMat)
  if (countsFromAbundance == "lengthScaledTPM") {
    newCounts <- abundanceMat * rowMeans(lengthMat)
  } else if (countsFromAbundance == "scaledTPM") {
    newCounts <- abundanceMat
  } else {
    stop("expecting 'lengthScaledTPM' or 'scaledTPM'")
  }
  newSum <- colSums(newCounts)
  countsMat <- t(t(newCounts) * (countsSum/newSum))
  countsMat
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

# code contributed from Andrew Morgan
read_kallisto_h5 <- function(fpath, ...) {
  if (!requireNamespace("rhdf5", quietly=TRUE)) {
    stop("reading kallisto results from hdf5 files requires Bioconductor package `rhdf5`")
  }
  counts <- rhdf5::h5read(fpath, "est_counts")
  ids <- rhdf5::h5read(fpath, "aux/ids")
  efflens <- rhdf5::h5read(fpath, "aux/eff_lengths")

  stopifnot(length(counts) == length(ids)) 
  stopifnot(length(efflens) == length(ids))

  result <- data.frame(target_id = ids,
                       eff_length = efflens,
                       est_counts = counts,
                       stringsAsFactors = FALSE)
  normfac <- with(result, (1e6)/sum(est_counts/eff_length))
  result$tpm <- with(result, normfac*(est_counts/eff_length))
  return(result)
}

summarizeFail <- function() {
  stop("

  tximport failed at summarizing to the gene-level.
  Please see 'Solutions' in the Details section of the man page: ?tximport

")
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

