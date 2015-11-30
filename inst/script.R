# read in some kallisto example files
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

library(devtools)
load_all()

txi <- tximport(files, type="kallisto", gene2tx=gene2tx)
head(txi$counts)

# generate counts from abundances and average transcript length over samples
# this way the counts are not correlated with length, and so the length matrix
# does not need to be provided as offset
txi2 <- tximport(files, type="kallisto", gene2tx=gene2tx, countsFromAbundance=TRUE)
head(txi2$counts)

#######################
# salmon
files <- file.path(dir,"salmon",samples$run,"quant.sf")
file.exists(files)

txi3 <- tximport(files, type="salmon", gene2tx=gene2tx)
head(txi3$counts)

#######################
# RSEM
files <- file.path(dir,"rsem",samples$run,paste0(samples$run,".genes.results"))
file.exists(files)

txi4 <- tximport(files, type="rsem")
head(txi4$counts)
