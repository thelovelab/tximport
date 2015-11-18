# put example code here
library(tximportData)
dir <- system.file("extdata",package="tximportData")
list.files(dir)
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
files <- file.path(dir,"kallisto",samples$run,"abundance.tsv")
file.exists(files)

gene2tx <- read.csv(file.path(dir, "gene2tx.csv"))

# these need gene IDs
# first make a table with columns: geneID, txID
lens <- readLengths(files, level="tx", gene2tx=gene2tx)
