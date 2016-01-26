# tximport 0.0.19

* Added `summarizeToGene` which breaks out the gene-level
  summary step, so it can be run by users on lists of
  transcript-level matrices produced by `tximport` with
  `txOut=TRUE`.

# tximport 0.0.18

* Changed argument `gene2tx` to `tx2gene`.
  This order is more intuitive: linking transcripts to genes,
  and matches the `geneMap` argument of Salmon and Sailfish.
