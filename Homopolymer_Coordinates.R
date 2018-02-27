###### Script to find homopolymer (poly-A/T/C/G) stretches in a given BSgenome with of user-defined minimum length.
###### Written by Alexander Toenges (a.toenges@uni-muenster.de)

###################################################################################################################
###################################################################################################################
require(BSgenome.Hsapiens.UCSC.hg38)
require(Biostrings)
require(GenomicRanges)
## Avoid floating point numbers for larger genomic coordinates:
options(scipen=999)
###################################################################################################################
###################################################################################################################

###################################################################################################################
###################################################################################################################

## Helper function for pattern matching:
## Function searches a given chromosome (chr) of the BSgenome (tmp.genome) for pattern matches (given.seq)
## Example: seq.check("AAAAAA", "chr1", BSgenome.Hsapiens.UCSC.hg38)
seq.check <- function(given.seq, chr, tmp.genome){
  
  ## Subset BSgenome to chromosome of interest:
  current.chr <- tmp.genome[[which(tmp.genome@seqinfo@seqnames == chr)]]
  
  ## Get all (also redundant) coordinates of the exact pattern match:
  tmp.match <- matchPattern(given.seq, getSeq(tmp.genome, chr))
  
  ## Transform to non-redundant granges:
  return( GenomicRanges::reduce( GRanges(seqnames = chr, ranges = ranges(tmp.match))) )
 
}

## Function that uses seq.check() to output a granges with coordinates of all defined homopolymers in the genome:
find.Homopolymers <- function(Nucleotide, NLength, Query.Genome){
  return(
    do.call("c", mclapply(paste("chr", c(seq(1,22), "X", "Y"), sep=""), function(x) 
      seq.check(
        paste(replicate(NLength, Nucleotide), collapse = ""), x, Query.Genome), mc.cores=8)
    )
  )
}

###################################################################################################################
###################################################################################################################

## Here an example for all homopolymers of length 6 or larger:
Homopolymer6_hg38.gr  <- do.call("c", list(find.Homopolymers("A", 6, BSgenome.Hsapiens.UCSC.hg38),
                                           find.Homopolymers("T", 6, BSgenome.Hsapiens.UCSC.hg38),
                                           find.Homopolymers("C", 6, BSgenome.Hsapiens.UCSC.hg38),
                                           find.Homopolymers("G", 6, BSgenome.Hsapiens.UCSC.hg38)
                                          )
                                 )

# Write 1-based granges to 0-based BED file on disk:
write.table(
  data.frame(
    seqnames(homoPolymerHexa.gr),
    start(homoPolymerHexa.gr)-1,
    end(homoPolymerHexa.gr)
  ), quote = F, col.names = F, row.names = F, sep="\t", file="Homopolymers_6AndMore_hg38.bed")

