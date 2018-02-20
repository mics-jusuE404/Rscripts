#### Script takes a BAM file which was loaded into R with GenomicALignments and writes it to bedGraph format on disk.
#### Optional RPM normalization and GZIP compression are possible.
#### Written by ALexander Toenges (a.toenges@uni-muenster.de), 2016

#############################################################################################################################
#############################################################################################################################

require(GenomicAlignments)

## Avoid floating point numbers for large chromosome coordinates
options(scipen=999)

#############################################################################################################################
#############################################################################################################################

## Tiny function:
GRanges2bedGraph <- function(GR, bedGraph, RPM, Gzip){
  
  totalReads   <- length(GR)
  tmp.Coverage <- as(coverage(GR), "GRanges")
  
  if(RPM == "y"){
    tmp.Coverage$score <- round((10^6 * tmp.Coverage$score * (1/totalReads)), digits = 6)
  }
  
  if(Gzip == "y") {
    write.table(
      data.frame(seqnames(tmp.Coverage), start(tmp.Coverage)-1, end(tmp.Coverage), tmp.Coverage$score),
      sep="\t", quote = F, col.names = F, row.names = F, file = gzfile(bedGraph)
    )
  } else {
      write.table(
        data.frame(seqnames(tmp.Coverage), start(tmp.Coverage)-1, end(tmp.Coverage), tmp.Coverage$score),
        sep="\t", quote = F, col.names = F, row.names = F, file = bedGraph
      )
    }
}


## Example
GRanges2bedGraph(GR = data.granges, bedGraph = "~/file.bedGraph.gz", RPM = "y", Gzip = "y")


