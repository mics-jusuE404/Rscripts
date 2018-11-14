#!/home/a/a_toen03/anaconda3/bin/Rscript

##############################################################################################################################

## Given a set of peaks and ATAC-seq data from two different conditions,
## infer differentially accessable motifs using chromVAR:

##############################################################################################################################

## Packages:
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(chromVARmotifs)
library(BSgenome.Mmusculus.UCSC.mm10)
library(data.table)
require(JASPAR2018)

set.seed(2018)
##############################################################################################################################

## Register cores for parallel processing:
register(MulticoreParam(16))

##############################################################################################################################

## Function to load BED files as GRanges():
Bed2GR <- function(BED){
  tmp.bed <- fread(BED, header = F, data.table = F)
  tmp_granges <- GRanges(seqnames=as.character(tmp.bed[,1]),
                         ranges=IRanges(
                           start = tmp.bed[,2]+1, end = tmp.bed[,3])
  ); return(tmp_granges)
}

##############################################################################################################################

## function (borrowed from chromVAR/TFBStools) to access JASPAR motifs (here the CORE vertebrate collection):
jaspar <- function (collection = "CORE", ...) 
{
  opts <- list()
  opts["tax_group"] <- "vertebrates"
  opts["collection"] <- collection
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out)))) 
    names(out) <- paste(names(out), TFBSTools::name(out), 
                        sep = "_")
  return(out)
}
motifs_JASPAR2018 <- jaspar()

##############################################################################################################################

## Step 1 -- Read in a list of reference ATAC-seq regions.
peaks_200bp <- readNarrowpeaks(filename = "merged_PeakSet.narrowPeak", 
                                         width = 200, 
                                         non_overlapping = T)

## Step 2 -- Get counts for a list of BAM files:
bam.list <- c("bam1.bam", "bam2.bam", "bam3.bam", "bam4.bam")
conditions <- c("wt", "wt", "treat1", "treat1")

peaks_counts          <- getCounts(alignment_files = c(bam.list), 
                                   peaks           = peaks_untreated_200bp, 
                                   paired          = TRUE, 
                                   by_rg           = FALSE, 
                                   format          = "bam", 
                                   colData         = DataFrame(condition = c("conditions")))

##############################################################################################################################

## Step 2 -- Define a function that runs the standard workflow, which is:
## - preprocessing of the count matrix (adding GC bias, filtering for depth, all default options),
## - matching motifs to the regions
## - computeDeviations and computeVariability (output is Outname_deviation/variability)
## - the direction of the variability output is determined by checking if the deviation score 
## --  if smaller or larger than zero upon the second condition.
## --  if smaller, than the variability score gets a *(-1), if larger, it stays positive

Bam2Deviation <- function(FragmentCounts, Peaks, Outname, Genome){
  
  ## this is pretty much copied from the chromVAR manual page:
  fragment_counts <- addGCBias(FragmentCounts, genome = Genome)
  filtered_counts <- filterSamples(fragment_counts, shiny = F)
  filtered_counts <- filterPeaks(fragment_counts)
  
  ## match motifs to the input regions:
  motif_ix <- matchMotifs(motifs, filtered_counts,
                                 genome = Genome)
  
  ## get deviation score:
  tmp.dev <- computeDeviations(object = filtered_counts, annotations = motif_ix)
  
  assign(paste(Outname, "_deviation", sep=""),
         tmp.dev,
         envir = .GlobalEnv
  )
  
  tmp.vari <- computeVariability(object = tmp.dev)
  
  ## if the deviation is negative, so less accessable in condition1, put a minus in front of the variablility score:
  tmp.z <- assays(tmp.dev)$z
  tmp.means <- sapply(unique(tmp.dev$celltype), function(x) rowMeans(tmp.z[,grep(x, tmp.dev$celltype)]))
  tmp.direction <- as.numeric(which( tmp.means[,1] > tmp.means[,2] ))
  
  tmp.vari$variability[tmp.direction] <- tmp.vari$variability[tmp.direction] * (-1)
  
  assign(paste(Outname, "_variability", sep=""),
         tmp.vari,
         envir = .GlobalEnv
  )
}

## Step-3 -- Run the function:
Bam2Deviation(FragmentCounts = peaks_counts, Peaks = peaks_200bp, Genome = BSgenome.Mmusculus.UCSC.mm10, Outname = "FOO")

## Write results to disk:
write.table(FOO_variability, sep="\t", quote = F, col.names = T, row.names = F, file="FOO_variability.tsv")

##############################################################################################################################
