#### Extract DNA sequence from a given BSgenome by coordinates and copy to clipboard
#### Tested only on macOS 10.12.5
#### Written by Alexander Toenges, 02/2018
#### last mod: 17.05.2018
#### BSgenome takes 1-based coordinates
#### ... and the script does not protects you from overloading your clipboard when querying a lot of basepairs, so be careful =)

##################################################################################################################################
##################################################################################################################################

require(Biostrings)
require(clipr)

##################################################################################################################################
##################################################################################################################################

#### Note: Coordinates must be supplied as "chr:start-end"

GetSeq2Clipboard <- function(BSgenome, Coords, Strand = "+"){
  
  if (Strand != "+" && Strand != "-") stop("Strand must either by + or -")
  chr <- strsplit(Coords, split=":")[[1]][1]
  tmp <- strsplit(
    strsplit(Coords, split=":")[[1]][2], split="-")[[1]]
  St <- tmp[1]
  En <- tmp[2]
  
  ## Save sequence as variable
  tmp.DNAstring <- toString(getSeq(BSgenome, chr, as.numeric(St), as.numeric(En)))
  
  ## If minus strand is specified, use reverseComplement:
  if (Strand == "-"){
    tmp.DNAstring <- as.character(reverseComplement(DNAString(tmp.DNAstring)))
  } 
  ## and spill to clipboard:
  clipr::write_clip(tmp.DNAstring)
}

## Example:
require(BSgenome.Hsapiens.UCSC.hg38)
GetSeq2Clipboard(BSgenome.Hsapiens.UCSC.hg38, "chr3:94646203-94646217", Strand = "-")
