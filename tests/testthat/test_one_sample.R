dir <- system.file("extdata", package="tximportData")
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
files <- file.path(dir,"salmon", samples$run, "quant.sf")
names(files) <- paste0("sample",1:6)
gene2tx <- read.csv(file.path(dir, "gene2tx.csv"))

txi <- tximport(files[1], type="salmon", gene2tx=gene2tx, reader=read_tsv)
expect_true(ncol(txi$counts) == 1)
