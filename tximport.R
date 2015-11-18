# put functions here
readLengths <- function(files,
                        level=c("tx","gene"),
                        geneIdCol="gene_id",
                        lengthCol="length",
                        abundanceCol="FPKM",
                        importer=read.table, ...) {
  dataList <- list()
  level <- match.arg(level, c("tx","gene"))
  # a weighted average of transcript length, weighting by transcript abundance
  avgTxLength <- function(x) {
    sum(x[[lengthCol]] * x[[abundanceCol]])/sum(x[[abundanceCol]])
  }
  for (i in seq_along(files)) {
    raw <- as.data.frame(importer(files[i], ...))
    if (level == "tx") {
      stopifnot(all(c(geneIdCol, lengthCol, abundanceCol) %in% names(raw)))
      raw[[geneIdCol]] <- factor(raw[[geneIdCol]], unique(raw[[geneIdCol]]))
      res <- do.call(c, as.list(by(raw, raw[[geneIdCol]], avgTxLength, simplify=FALSE)))      
      dataList[[i]] <- data.frame(avgLength=res, row.names=names(res))
    } else if (level == "gene") {
      stopifnot(all(c(geneIdCol, lengthCol) %in% names(raw)))
      dataList[[i]] <- data.frame(avgLength=raw[[lengthCol]], row.names=raw[[geneIdCol]])
    }
  }
  data <- do.call(cbind, dataList)
  as.matrix(data)
}
