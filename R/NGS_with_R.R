
install_package <- function() {
  source("http://bioconductor.org/biocLite.R")
  source("http://bioconductor.org/scratch-repos/pkgInstall.R")
  pkgInstall("SequenceAnalysis")
}

library(BSgenome.Hsapiens.UCSC.hg19)
?BSgenome.Hsapiens.UCSC.hg19

genome <- BSgenome.Hsapiens.UCSC.hg19
seqlengths(genome)
plot(seqlengths(genome))
seqnames(genome)
genome$chr1  # same as genome[["chr1"]]

for (seqname in seqnames(genome)[1:2]) {
  cat("Checking sequence", seqname, "... ", length(genome[[seqname]]))
  #seq <- genome[[seqname]]
  #checkOnlyNsInGaps(seq)
  cat("\n")
}

length(names(genome$upstream1000))
names(genome$upstream1000)[1:4]

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)


getSNPcount()
plot(getSNPcount(),seqlengths(genome)[1:25])

