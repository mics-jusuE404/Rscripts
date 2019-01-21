## Simply function to check a given DNA sequence (e.g. a primer sequence) for perfect match 
## against both strands of a chromosome from a BSgenome, returning a data frame with the 
## position and strand (or an empty one if no match):

CheckSequence <- function(QUERY, CHR, GENOME){
  
  library(Biostrings)
  
  ## select chromosome of BSgenome:
  toScan <- GENOME[[which(GENOME@seqinfo@seqnames == CHR)]]
  
  ## match sequence and reverse complement:
  matched_fwd <- matchPattern(QUERY, toScan)
  matched_rev <- matchPattern(as.character(reverseComplement(DNAString(QUERY))), toScan)
 
   ## if no result:
  if (length(matched_fwd) == 0 & length(matched_rev) == 0){
    df.out <- data.frame(chr="NA", start="NA", end="NA", sequence=QUERY, strand="NA")
  }
  
  ## if result:
  df.out <- data.frame(chr=NULL, start=NULL, end=NULL, sequence=NULL, strand=NULL)
  if (length(matched_fwd) > 0) df.out <- rbind(df.out, data.frame(chr=CHR, start=start(matched_fwd), end=end(matched_fwd), sequence=QUERY, strand=as.character("+")))
  if (length(matched_rev) > 0) df.out <- rbind(df.out, data.frame(chr=CHR, start=start(matched_rev), end=end(matched_rev), sequence=QUERY, strand="-"))
  return(df.out)
}

## Example for a single chromosome:
library(BSgenome.Hsapiens.UCSC.hg38)
CheckSequence(QUERY = "CAACAAGGTGCCAAGTCTTTT", CHR = "chr11", GENOME = BSgenome.Hsapiens.UCSC.hg38)

## and for all chromosomes (takes like 10 seconds or so):
do.call(rbind, 
        mclapply(paste("chr", c(seq(1,22), "X", "Y"), sep=""), function(x) {
          return(CheckSequence(QUERY = "CAACAAGGTGCCAAGTCTTTT", CHR = x, GENOME = BSgenome.Hsapiens.UCSC.hg38))
          }, mc.cores=16
        )
)

