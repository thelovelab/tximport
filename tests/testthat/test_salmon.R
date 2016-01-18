dir <- system.file("extdata", package="tximportData")
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
files <- file.path(dir,"salmon", samples$run, "quant.sf")
names(files) <- paste0("sample",1:6)
tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))

txi <- tximport(files, type="salmon", tx2gene=tx2gene, reader=read_tsv)
expect_true(ncol(txi$counts) == length(files))

# also test txOut here
txi.txout <- tximport(files, type="salmon", txOut=TRUE, reader=read_tsv)
expect_true(ncol(txi.txout$counts) == length(files))

# test reading in slow way
txi <- tximport(files[1:2], type="salmon", tx2gene=tx2gene)
expect_true(ncol(txi$counts) == 2)
