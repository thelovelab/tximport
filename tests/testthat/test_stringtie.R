context("stringtie")
test_that("import stringtie works", {

  # these files created with the command:
  # stringtie -eB -G chess1.0.gff <source_file.bam> 

  makeData <- function(n) {
    data.frame(t_id = 1:n,
               chr = rep("chr1",n),
               strand = rep("+",n),
               start = 1:n * 1e4,
               end = 1:n * 1e4 + 1000,
               t_name = 1:n,
               num_exons = rep(10,n),
               length = 1:n * 1e3,
               gene_id = rep(1:10,each=n/10),
               gene_name = rep(letters[1:10],each=n/10),
               cov = rpois(n, 100),
               FPKM = rnorm(n, 100, 10))
  }
  n <- 30
  A <- makeData(n)
  B <- makeData(n)
  C <- makeData(n)
  files <- c(A="A", B="B", C="C")
  importer <- function(x) get(x)
  tx2gene <- A[,c("t_name","gene_name")]
  txi <- tximport(files, type="stringtie", tx2gene=tx2gene,
                  importer=importer, existenceOptional=TRUE)

  skip_on_os("windows")
  
  expect_true(txi$counts[1,1] == sum(A$cov[1:3] * A$length[1:3] / 75))

})
