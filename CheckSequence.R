## Simply function to checks a BSgenome for occurrence of a given sequence, returning a reduced GRanges:

###### Script to find a specific pattern in a BSgenome returning a GRanges object with the coordinates:

require(Biostrings)
require(GenomicRanges)
require(BSgenome.Hsapiens.UCSC.hg38)
options(scipen=999)

FindPolyX <- function(QUERY, BSGENOME, CORES=16){

  all.matches <- mclapply(1:length(seqnames(BSGENOME)), function(x) {

                    tmp.match <- matchPattern(QUERY, BSGENOME[[x]])

                    if (length(tmp.match) > 0) {
                      tmp.gr    <- GRanges(seqnames = seqnames(BSGENOME)[x],
                                                      ranges = ranges(tmp.match))
                      return(GenomicRanges::reduce(tmp.gr))
                    }

                  }, mc.cores = CORES)
  return(suppressWarnings( do.call("c", all.matches) ) )
}

FindPolyX(QUERY = "AAAAAAAAAAAA", BSGENOME = BSgenome.Hsapiens.UCSC.hg38)
