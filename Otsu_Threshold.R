#### Take a list of integers and find a threshold that splits these integers into two groups,
#### aiming to minimize intragroup variance but maximize intergroup variance.
#### Based on Otsu's threshold and code pretty much copied from <https://en.wikipedia.org/wiki/Otsu%27s_method>
#### Modifications are the log10 transformation to make it suitable for larger numbers as common in NGS.
#### Might be suitable if NGS data show strong bimodal distribution and any thresholding is desired (beyond by-eye methods)

Otsu_NGS <- function(COUNTS){
  tmp.table <- table(log10(COUNTS + 1))
  tmp.tabname <- as.numeric(attr(tmp.table,"dimnames")[[1]])
  tmp.occur <- as.vector(tmp.table)
  tmp.density <- tmp.occur/sum(tmp.occur) #probs (tabreadcount/total)
  
  # Set everything to zero:
  maxi <- 0; sum1 <- 0; sum1 <- sum(tmp.tabname * tmp.density); sumB <- 0
  wB <- 0;  wF <- 0
  
  ## Iterate because iterations are fun!
  for (q in 1: length(tmp.tabname)){
    wB <- wB + tmp.density[q]
    if (wB == 0) message("wB ist 0") 
    wF <- 1 - wB
    if (wF == 0) break
    sumB <- sumB + tmp.tabname[q] * tmp.density[q]
    mB <- sumB/wB
    mF <- (sum1-sumB)/ wF
    interVar <- wB*wF*(mB-mF)*(mB-mF)
    if (interVar > maxi) {
      maxi.level <- tmp.tabname[q] #maxi.level is the greatest value that still belongs to the background/lower class
      maxi <- interVar
    }
  }
  return (round(10^maxi.level) - 1)
}
