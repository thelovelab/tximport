dir <- system.file("extdata", package="tximportData")
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
files <- file.path(dir,"kallisto", samples$run, "abundance.tsv")
names(files) <- paste0("sample",1:6)
gene2tx <- read.csv(file.path(dir, "gene2tx.csv"))

txi.cfa <- tximport(files, type="kallisto", gene2tx=gene2tx, 
                    countsFromAbundance="scaledTPM", reader=read_tsv)
txi.cfa.len <- tximport(files, type="kallisto", gene2tx=gene2tx, 
                        countsFromAbundance="lengthScaledTPM", reader=read_tsv)

expect_true(ncol(txi.cfa$counts) == length(files))
expect_true(ncol(txi.cfa.len$counts) == length(files))
