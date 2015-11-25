# read in some example files
library(tximportData)
dir <- system.file("extdata",package="tximportData")
list.files(dir)
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
files <- file.path(dir,"kallisto",samples$run,"abundance.tsv")
file.exists(files)

# transcripts need to be associated with gene IDs
# if that info is in the files, skip this step
# if not, first make a data.frame with two columns: gene ID, transcript ID
gene2tx <- read.csv(file.path(dir, "gene2tx.csv"))

system.time({ res <- tximport(files, level="tx", gene2tx=gene2tx) })
