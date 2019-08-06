## Check BSgenome for perfect matches of a nucleotide string,
## outputs reduced GRanges object:

CheckSequence <- function(QUERY, BSGENOME, CORES=16){

  require(Biostrings)
  require(GenomicRanges)
  options(scipen=999)
  
  ## remove whitespaces if input contains some:
  QUERY < -DNAString( gsub(" ", "", as.character(QUERY)) )
  
  all.matches <- mclapply(1:length(seqnames(BSGENOME)), function(x) {

                    tmp.match <- matchPattern(QUERY, BSGENOME[[x]])

                    if (length(tmp.match) > 0) {
                      tmp.gr    <- GRanges(seqnames = seqnames(BSGENOME)[x],
                                                      ranges = ranges(tmp.match))
                      return(GenomicRanges::reduce(tmp.gr))
                    }

                  }, mc.cores = CORES)
  return(suppressWarnings( do.call("c", unlist(all.matches) ) ) )
}

library(BSgenome.Hsapiens.UCSC.hg38)
CheckSequence(QUERY = "AAAAAAAAAAAA", BSGENOME = BSgenome.Hsapiens.UCSC.hg38)
