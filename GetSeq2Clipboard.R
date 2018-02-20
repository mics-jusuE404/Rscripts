#### Extract DNA sequence from a given BSgenome by coordinates and copy to clipboard
#### Tested only on macOS 10.12.5
#### Written by Alexander Toenges, 02/2018
#### BSgenome takes 1-based coordinates
#### ... and the script does not protects you from overloading your clipboard when querying like 10^6 basepairs, so be careful =)

##################################################################################################################################
##################################################################################################################################

require(BSgenome.Hsapiens.UCSC.hg38)
require(Biostrings)
require(clipr)

##################################################################################################################################
##################################################################################################################################

#### Note: Coordinates must be supplied as "chr:start-end"

GetSeq2Clipboard <- function(BSgenome, Coords){
  
  chr <- strsplit(Coords, split=":")[[1]][1]
  tmp <- strsplit(
          strsplit(Coords, split=":")[[1]][2], split="-")[[1]]
  St <- tmp[1]
  En <- tmp[2]
  clipr::write_clip(toString(getSeq(BSgenome, chr, as.numeric(St), as.numeric(En))))
}

## Example:
GetSeq2Clipboard(BSgenome.Hsapiens.UCSC.hg38, "chr1:444444-444460") # AATCCAGGAGCTGGTTT
