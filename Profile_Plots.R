#### Script to make profile plots (= read intensities of e.g. ChIP-seq) over a used-supplied list of regions.
#### Written by Alexander Toenges (a.toenges@uni-muenster.de), 2017
#### There is still plenty of room to automate tasks in this script:

#############################################################################################################################
#############################################################################################################################

#### Script assumes that a normalized BigWig file of the ChIP-seq or whatever-experiment has been created, e.g. with DEEPtools,
#### and also needs a set of regions as reference in BED format, both from disk.

#############################################################################################################################
#############################################################################################################################

require(genomation)
require(data.table)
require(GenomicRanges)

## Load reference regions:
regions_foreground.bed <- fread("~/IMTB/DLBCL/TAC-STARRseq/pSTARRseq_human/180115_enhancer.bed", header = F, data.table = F)
regions_background.bed <- fread("~/IMTB/DLBCL/TAC-STARRseq/pSTARRseq_human/180115_NOenhancer.bed", header = F, data.table = F)

## Write BED files as GRanges and define a window around them:
Bed2GRanges_WithWindow <- function(dataframe, winsize){
  if( class(dataframe) != "data.frame" ) stop("Please load regions as data-frame. Check if data.table=F")
  
  tmp.gr <- GRanges(seqnames = dataframe[,1],
                    ranges = IRanges(start = dataframe[,2]+1,
                                     end   = dataframe[,3])
  )
  tmp.centers <- start(tmp.gr) + ceiling(0.5*width(tmp.gr))
  # Windows around the centers:
  tmp.windows <- GRanges(seqnames = seqnames(tmp.gr),
                         ranges = IRanges(start = tmp.centers - winsize +1 ,
                                          end   = tmp.centers + winsize)
  )
  return (tmp.windows)
}

## The first GRanges is the foreground, so regions of interest:
foreground.window <- Bed2GRanges_WithWindow(dataframe = regions_foreground.bed, winsize = 10000)

## The second one is a set of background regions, e.g. from a genome-wide reporter assay, regions that turned out to be 
## negative in the assay, where foreground = regions that turned out positive.
## the sample() takes an equal number of backgrounds to match the number of foregrounds:
background.window <- Bed2GRanges_WithWindow(
  dataframe = regions_background.bed[sample(seq(1:nrow(background.window)), size = nrow(foreground.window)),],
  winsize = 10000)
  
#############################################################################################################################
#############################################################################################################################

#### Calculate a genomation::ScoreMatrix().
#### Assumes that bigWigs are already normalized for read depth:
Score.foreground <- ScoreMatrix(target = "/data.bigwig_RPM.bigwig", type = "bigWig",
                               windows = foreground.window, strand.aware = F, rpm = F, unique = F)
                               
Score.background <- ScoreMatrix(target = "/data.bigwig_RPM.bigwig", type = "bigWig",
                                 windows = background.window, strand.aware = F, rpm = F, unique = F)
                                 
#############################################################################################################################
#############################################################################################################################

#### Plot the data:
plot.Profile <- function(foreground.data, background.data, main){
  
  ## Check which score matrix has the highest value so that this is plotted first to set the upper ylim:
  if (max(colMeans(foreground.data)) > max(colMeans(background.data))){
    Larger <- foreground.data
    Smaller <- background.data
  } else {
    Larger <- background.data
    Smaller <- foreground.data
  }
  
  # Plot the data
  plot(
    colMeans(Larger), type="l", bty="n", ylab="Average RPM", xlab="", col="firebrick", lwd=1, xaxt="n", main=main,
       ylim=c(min(colMeans(Smaller)) * 1.1, max(colMeans(Larger)*1.1))
  )
  lines(colMeans(Smaller), lwd=1, col="darkblue")  
  
  wsize <- dim(foreground.data)[2] / 2
  lower <- paste("-", 0.001*dim(foreground.data)[2]/2, "kb", sep="")
  upper <- paste("+", 0.001*dim(foreground.data)[2]/2, "kb", sep="")
  
  axis(side = 1, at = c(1, ncol(foreground.data)/2, ncol(foreground.data)), labels=c(lower, "0kb", upper)) 
   
  legend("topright", legend=c("Enhancers", "Background"), bty="n", lty=1, lwd=1, col=c("firebrick", "darkblue"))
}

## Plot:
plot.Profile(foreground.data = Score.foreground, background.data = Score.background, main = "Plot Title")

#### Should then look e.g. like this: https://uni-muenster.sciebo.de/s/bVecCBzyMGoShjm
                                 
