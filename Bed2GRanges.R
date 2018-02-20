#### Read in a BED file and write as GRanges.
#### Written by Alexander Toenges (a.toenges@uni-muenster.de)
#### ChromSizes is a file with $1 = chromosome names and $2 the length of it, tab-delim
#### => can be downloaded from UCSC together with the reference genome fasta files

########################################################################################################################
########################################################################################################################

Bed2GRanges <- function(BED, ChromSizes, Stranded){
  
  require(GenomicRanges)
  require(data.table)
  
  tmp.bed <- fread(file = BED, data.table = F)
  tmp.chromSize <- read.table(file = ChromSizes, sep="\t", header = FALSE)
  
  ## GRanges constructor: Note that GRanges is 1-based and BED is 0-based, therefore tmp.bed[,2] + 1:
  tmp_granges <- GRanges(seqnames=as.character(tmp.bed[,1]),
                         ranges=IRanges(
                           start = tmp.bed[,2]+1, end = tmp.bed[,3])
                         )
  ## Add optional name information if $4 in BED exists
  if (ncol(tmp.bed) > 3) names(tmp_granges) <- tmp.bed[,4]
  
  ## Add optional strand information:
  if (Stranded == "y" && ncol(tmp.bed) >= 6) strand(tmp_granges) <- tmp.bed[,6]
  
  ## Add seqlengths (= length of the chromosomes) according to chromSizes file:
  foo <- c()
  for (i in seqnames(tmp_granges)@values) {
    foo <- c(foo,which(i == tmp.chromSize[,1]))
  }
  seqlengths(tmp_granges) <- tmp.chromSize[,2][foo]
  
  ## Trim GRanges in case some entries are out-of-bounds according to chromSizes/seqlengths:
  return(trim(tmp_granges))
}

## Example:
Bed2GRanges(BED = "~/test.bed", ChromSizes = "/Volumes/Rumpelkammer/Genomes/hg19/hg19_chromSizes.txt", Stranded = "y")

