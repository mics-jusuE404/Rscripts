## Take a list of integers and apply Otsu's method to find a threshold that 
## maximizes inter-group variance but minimizes intra-group variance.

Otsu_NGS <- function(COUNTS, log10 = TRUE){
  
  if (log10) tmp.table    <- table(log10(COUNTS + 1))
  if (!log10) tmp.table   <- table(COUNTS + 1)
  
  tmp.tabname <- as.numeric(attr(tmp.table,"dimnames")[[1]])
  tmp.occur   <- as.vector(tmp.table)
  tmp.density <- tmp.occur/sum(tmp.occur) #probs (tabreadcount/total)
  
  # Set everything to zero:
  maxi <- 0; sum1 <- 0; sum1 <- sum(tmp.tabname * tmp.density); sumB <- 0
  wB   <- 0;  wF  <- 0
  
  ## Iterate because iterations are fun!
  for (q in 1: length(tmp.tabname)){
    wB <- wB + tmp.density[q]
    if (wB == 0) message("wB ist 0") 
    wF <- 1 - wB
    if (wF == 0) break
    sumB <- sumB + tmp.tabname[q] * tmp.density[q]
    mB   <- sumB/wB
    mF   <- (sum1-sumB)/ wF
    interVar <- wB*wF*(mB-mF)*(mB-mF)
    if (interVar > maxi) {
      maxi.level <- tmp.tabname[q] # maxi.level is the largest value that still belongs to the lower group
      maxi <- interVar
    }
  }
  if (log10)  return (round(10^maxi.level) - 1)
  if (!log10) return (maxi.level)
}
