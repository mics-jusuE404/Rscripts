###### Script to find homopolymer (poly-A/T/C/G) stretches in a given BSgenome with of user-defined minimum length.
###### Written by Alexander Toenges (a.toenges@uni-muenster.de)

###################################################################################################################

require(Biostrings)
require(GenomicRanges)

###################################################################################################################

## Avoid floating point numbers for larger genomic coordinates:
options(scipen=999)

###################################################################################################################

## Helper function for pattern matching:
## Function searches a given chromosome (chr) of the BSgenome (tmp.genome) for pattern matches (given.seq)
## Example: seq.check("AAAAAA", "chr1", BSgenome.Hsapiens.UCSC.hg38)

seq.check <- function(given.seq, chr, tmp.genome){
  
  ## Subset BSgenome to chromosome of interest:
  current.chr <- tmp.genome[[which(tmp.genome@seqinfo@seqnames == chr)]]
  
  ## Get all (also redundant) coordinates of the exact pattern match:
  tmp.match <- matchPattern(given.seq, getSeq(tmp.genome, chr))
  
  ## If no matches are found, output empty GRanges, else output GRanges with unique coordinates:
  if (length(tmp.match@ranges) == 0){
    return(
      GRanges()
    )
  } else{ 
      return( 
        suppressWarnings(GenomicRanges::reduce(GRanges(seqnames = chr, ranges = ranges(tmp.match))))
      )
    }
}

## Function that uses seq.check() to output a granges with coordinates of all defined homopolymers in the genome:
find.Homopolymers <- function(Nucleotide, PolyXLength, Query.Genome, Cores = 1){
  
  ## Lapply-based function to scan all chromosomes of the BSgenome for the intended Nucleotide of length PolyXLength
  ## Return as GRanges
  return(
    suppressWarnings(
      do.call("c", 
        mclapply(as.character(Query.Genome@seqinfo@seqnames), function(x) seq.check(
                      paste(replicate(PolyXLength, Nucleotide), collapse = ""), x, Query.Genome), 
                      mc.cores = Cores)
      )
    )
  )
}

###################################################################################################################

## Example:
require(BSgenome.Ecoli.NCBI.20080805)

## 1: Find polyA's of at least 5bp length:
PolyA5_EColi <- find.Homopolymers(Nucleotide = "A", PolyXLength = 5, Query.Genome = BSgenome.Ecoli.NCBI.20080805)

## 2: Find all (polyATCG) of at least 5bp length:
PolyAll5_EColi  <- do.call("c", list(find.Homopolymers("A", 5, BSgenome.Ecoli.NCBI.20080805),
                                     find.Homopolymers("T", 5, BSgenome.Ecoli.NCBI.20080805),
                                     find.Homopolymers("C", 5, BSgenome.Ecoli.NCBI.20080805),
                                     find.Homopolymers("G", 5, BSgenome.Ecoli.NCBI.20080805)
                                )
                    )

