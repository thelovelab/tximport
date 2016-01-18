dir <- system.file("extdata", package="tximportData")
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
files <- file.path(dir,"salmon", samples$run, "quant.sf")
names(files) <- paste0("sample",1:6)
tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))

txi <- tximport(files[1], type="salmon", tx2gene=tx2gene, reader=read_tsv)
expect_true(ncol(txi$counts) == 1)
