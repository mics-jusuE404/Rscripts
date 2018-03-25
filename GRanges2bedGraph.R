#### Script takes a BAM file which was loaded into R with GenomicALignments and writes it to bedGraph format on disk.
#### Optional RPM normalization and GZIP compression are possible.
#### Written by ALexander Toenges (a.toenges@uni-muenster.de), 2016

#############################################################################################################################
#############################################################################################################################

require(GenomicAlignments)

## Avoid floating point numbers for large chromosome coordinates which would corrupt the BED format
options(scipen=999)

#############################################################################################################################
#############################################################################################################################

## Wrapper:
GRanges2bedGraph <- function(GR, bedGraph, RPM = "n", Gzip = "n"){
  
  if (class(GR) != "GRanges")  stop("GR is not a GRanges object!")
  
  totalReads   <- length(GR)
  tmp.Coverage <- as(coverage(GR), "GRanges")
  
  if(RPM == "y"){
    tmp.Coverage$score <- round((10^6 * tmp.Coverage$score * (1/totalReads)), digits = 6)
  }
  
  if(Gzip == "y") {
    tmp.file <- gzfile(bedGraph)
  } else {
    tmp.file <- bedGraph
  }
  
  write.table(
    data.frame(seqnames(tmp.Coverage), start(tmp.Coverage)-1, end(tmp.Coverage), tmp.Coverage$score),
    sep="\t", quote = F, col.names = F, row.names = F, file = tmp.file
  )
}

## Example
GRanges2bedGraph(GR = data.granges, bedGraph = "~/file.bedGraph.gz", RPM = "y", Gzip = "y")


