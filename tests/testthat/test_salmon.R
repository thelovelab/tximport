dir <- system.file("extdata", package="tximportData")
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
files <- file.path(dir,"salmon", samples$run, "quant.sf")
names(files) <- paste0("sample",1:6)
gene2tx <- read.csv(file.path(dir, "gene2tx.csv"))

txi <- tximport(files, type="salmon", gene2tx=gene2tx, reader=read_tsv)
expect_true(ncol(txi$counts) == length(files))

# also test txOut here
txi.txout <- tximport(files, type="salmon", txOut=TRUE, reader=read_tsv)
expect_true(ncol(txi.txout$counts) == length(files))

# test reading in slow way
txi <- tximport(files[1:2], type="salmon", gene2tx=gene2tx)
expect_true(ncol(txi$counts) == 2)
