# put functions here
readLengths <- function(files,
                        level=c("tx","gene"),
                        geneIdCol="gene_id",
                        txIdCol="target_id",
                        lengthCol="eff_length",
                        abundanceCol="tpm",
                        gene2tx=NULL,
                        importer=function(x) read.table(x,header=TRUE),
                        ...) {
  dataList <- list()
  level <- match.arg(level, c("tx","gene"))
  # a weighted average of transcript length, weighting by transcript abundance
  avgTxLength <- function(x) {
    sum(x[[lengthCol]] * x[[abundanceCol]])/sum(x[[abundanceCol]])
  }
  for (i in seq_along(files)) {
    cat(i,"")
    raw <- as.data.frame(importer(files[i], ...))
    if (level == "tx") {
      # does the table contain gene association or was an external gene2tx table provided?
      if (is.null(gene2tx)) {
        stopifnot(all(c(geneIdCol, lengthCol, abundanceCol) %in% names(raw)))
      } else {
        stopifnot(all(c(lengthCol, abundanceCol) %in% names(raw)))
        # need that the transcripts are associated with genes
        raw <- raw[raw[[txIdCol]] %in% gene2tx$TXNAME, ]
        raw[[geneIdCol]] <- gene2tx$GENEID[match(raw[[txIdCol]], gene2tx$TXNAME)]
      }
      raw[[geneIdCol]] <- factor(raw[[geneIdCol]], unique(raw[[geneIdCol]]))
      res <- do.call(c, as.list(by(raw, raw[[geneIdCol]], avgTxLength, simplify=FALSE)))
      dataList[[i]] <- res
    } else if (level == "gene") {
      stopifnot(all(c(geneIdCol, lengthCol) %in% names(raw)))
      dataList[[i]] <- raw[[lengthCol]]
      names(dataList[[i]]) <- raw[[geneIdCol]]
    }
  }
  cat()
  as.matrix(do.call(cbind, dataList))
}
