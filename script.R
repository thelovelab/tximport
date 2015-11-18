# put example code here
library(tximportData)
dir <- system.file("extdata",package="tximportData")
list.files(dir)
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
files <- file.path(dir,"kallisto",samples$run,"abundance.tsv")
file.exists(files)

# these need gene IDs
# first make a table with columns: geneID, txID
lens <- readLengths(files, level="tx") 
