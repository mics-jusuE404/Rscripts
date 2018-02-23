#### A little dirty but functional script that reads the default output of annotatePeaks.pl from HOMER
#### and plots the number of peaks per genomic feature as a pie chart:
#### Author: Alexander Toenges (2016) a.toenges@uni-muenster.de

PieChart_annotatePeaksHomer <- function(AnnotFile, CustomTitle){
  
  require(data.table)
  
  ## Load homer output:
  AnnotFile <- fread(AnnotFile, data.table=FALSE, header=TRUE)
  ## Crude check if the first colname if OK
  if(grepl('PeakID', colnames(AnnotFile)[1]) != "TRUE") {
    stop("Input file does not appear to be output of annotatePeaks.pl")
  }
  
  ## By default homer calls everything promoter that is -1000 to +100 bp
  ## Capture the promoter entries and split them into -100/+100 by and -1kb to -0.1kb:
  
  file.prom <- AnnotFile[grepl("promoter", AnnotFile$Annotation),]
  file.rest <- AnnotFile[!grepl("promoter", AnnotFile$Annotation),]
  
  ## more = more than 0.1kb away from peak center, less = within
  tmp.more <- file.prom[file.prom$`Distance to TSS` > 100 | file.prom$`Distance to TSS` < -100,]
  tmp.more$Annotation <- c("TSS_-1000")
  tmp.less <- file.prom[file.prom$`Distance to TSS` <= 100 | file.prom$`Distance to TSS` >= -100,]
  tmp.less$Annotation <- c("TSS_-100/+100")
  
  ## bidn back to new df
  tmp.annot <- rbind(file.rest, tmp.more, tmp.less)
  
  a <- table(unlist(lapply(1:nrow(tmp.annot), function(x) strsplit(tmp.annot$Annotation[x], split=" ")[[1]][1])))
  
  names(a)[which(names(a) == "TSS_-100/+100")] <- c("TSS -100/+100bp")
  names(a)[which(names(a) == "TSS_-1000")] <- c("TSS -1000bp")
  
  a.names <- attr(a, "dimnames")[[1]]
  a.counts <- as.vector(a)
  
  ## rename 3' and 5' in 3'UTR and 5'UTR if these entries exist:
  if ( length(which(a.names == "3'")) > 0 ) {
    a.names[which(a.names == "3'")] <- "3'UTR"
  }
  if ( length(which(a.names == "5'")) > 0 ) {
    a.names[which(a.names == "5'")] <- "5'UTR"
  }
  
  ## Merge "non-coding" and "exon" into exon
  if ( ("non-coding" %in% a.names) == TRUE ){ #check if "non-coding" is present
    
    ## if exon is present, merge both
    if( ("exon" %in% a.names) == TRUE){
      merged.exonCount <- as.numeric(a.counts[which(a.names == "exon")]) + as.numeric(a.counts[which(a.names == "non-coding")])
      ## replace:
      a.counts[which(a.names == "exon")] <- merged.exonCount
    } else{ ## if exon is not present, simply rename "non-coding" to exon
      a.names[which(a.names == "non-coding")] <- "exon"
    }
  }
  if (("Intergenic" %in% a.names) == TRUE){
    a.names[which(a.names == "Intergenic")] <- "intergenic"
  }
  df.out <- data.frame(a.names, a.counts)
  df.out <- df.out[-which(df.out$a.names == "non-coding"),] 
  
  ## Get numbers per feature and plot as Pie:
  my.percent <- round(100*(as.numeric(df.out[,2]/sum(as.numeric(df.out[,2])), digits=1)))
  my.cols <- c("darkorchid4", "darkorange3", "black", "skyblue4", "goldenrod3", "red3", "darkolivegreen4", "cadetblue")
  my.label <- as.vector(df.out[,1])
  do.plot <- pie(as.numeric(df.out[,2]), labels = "", col = my.cols, lty = 1, lwd=1, radius = 0.56)
  title(CustomTitle, line= -3)
  legend(-0.58,-0.6, legend=sapply(my.label, as.expression), fill=my.cols, 
         bty="n", horiz = FALSE, y.intersp = 1.4, x.intersp = 0.5, text.font=1, cex=1, ncol = 2, text.width=0.8)
  
  return(c(df.out, do.plot))
}

## Example:
example.pie <- PieChart_HomerAnnotation(AnnotFile = "~/peaks_annotated.txt", CustomTitle = "ATAC-seq Peak Annotation")
## => Looks like this: https://uni-muenster.sciebo.de/s/tpMfJyDZgfypyZR

