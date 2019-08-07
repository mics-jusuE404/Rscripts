## Permutation input for InTAD:

## Hilariously long script that could probably be done with a few commands,
## but lets have some fun:
makeUniqueCombinations <- function(list.rna, list.atac){
  
  require(parallel)
  
  ## in case of trouble, stop without this ugly stop() error message,
  ## credit: https://stackoverflow.com/questions/14469522/stop-an-r-program-without-error
  stopQuietly <- function() {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }
  
  ## more than two cell types: not implemented:
  if (length( length(unique(sapply(strsplit(list.rna, split="_"), function(x)x[2]))) ) > 2){
    message("[Error] More than two cell types, currently not supported!")
    stopQuietly()
  }
    
  ## https://www.r-bloggers.com/outersect-the-opposite-of-rs-intersect-function/
  outersect <- function(x, y) {
    sort(c(setdiff(x, y),
           setdiff(y, x)))
  }
  
  if (0 %in% sapply(c(list.rna, list.atac), function(voltorb) nchar(voltorb))) {
    message("[Error] Check the input lists, seems that at least one contains empty strings or other odd entries")
    stopQuietly()
  }
  
  ## Get the unique cell types:
  unique.celltype.rna  <- unique(sapply(strsplit(list.rna, split="_"), function(x)x[2]))
  unique.celltype.atac <- unique(sapply(strsplit(list.atac, split="_"), function(x)x[2]))
  
  if ( length( outersect(unique.celltype.atac, unique.celltype.rna)) ){
    message("[Error] One list contains a cell type not present in the second.")
    stopQuietly()
  }
  
  combn.celltypes <- mclapply(unique.celltype.atac, function(G){
    
    tmp.list.rna  <- grep(G, list.rna, value = T)
    tmp.list.atac <- grep(G, list.atac, value = T)
    
    ## put some unique letters around the input names to ensure they are unique:
    tmp.list.rna <- paste0("RNASEQ", tmp.list.rna)
    tmp.list.atac <- paste0("ATACSEQ", tmp.list.atac)
    
    ## make sure elements are unique in each list:
    if ( length(tmp.list.rna) != length(unique(tmp.list.rna)) | 
         length(tmp.list.atac) != length(unique(tmp.list.atac)) ){
      message("[Error]: At least one of the two lists contains duplicated entries!")
      stopQuietly()
    }
    
    ## make all possible pairwise combinations
    ## https://stackoverflow.com/questions/13018173/elementwise-combination-of-two-lists-in-r
    combos.all <- outer(tmp.list.rna,tmp.list.atac, FUN=paste)
    
    ## separate by "__" as this is pretty sure a unique delimiter
    combos.all <- as.character(sapply(c(combos.all), function(x) gsub(" ", "__",x)))
    
    ## if less than three combinations (so only 1 element in either tmp.list.rna or tmp.list.atac),
    ## return data at this point as combos.all is already the full list of valid combinations:
    if (length(combos.all) < 3){
     
      combos.all <- gsub("RNASEQ|ATACSEQ", "", paste0("RNASEQ", combos.all))
      return(combos.all)
      
    } 
    
    ## => combos.all contains all possible combinations between the ATAC-seq and RNA-seq colnames
    ## per cell type.
    ## Now find all the allowed combinations within combos.all
    
    ## tmp.shortest is the number of elements we one per combination,
    ## so if cell type has 3 replicate, we want 3
    tmp.shortest <- min(sapply(c(list.atac, list.rna), function(x) {
      grep1 <- grep(G, list.atac)
      grep2 <- grep(G, list.rna)
      return(min(length(grep1), length(grep2)))
    }))
    
    ## make all possible combinations from combos.all no matter if valid
    combis <- combn(seq(1, length(combos.all)), tmp.shortest)
    
    ## now keep only those that are valid by checking that no sample of 
    ## either RNA or ATAC would occur > once in that combination
    tmp.ok <- sapply(1:dim(combis)[2], function(U){
      tmp.toExtract <- combis[,U]
      tmp.valid <- sort(combos.all[tmp.toExtract])
      if (sum(table(unlist(strsplit(tmp.valid, split="__"))) > 1) == 0) {
        return(paste(sort(combos.all[tmp.toExtract]), collapse="@@"))
      }
    })
    tmp.ok <- unlist(tmp.ok)
    tmp.ok <- gsub("RNASEQ|ATACSEQ", "", tmp.ok)
    return(sort(tmp.ok))
  })
  
  if (length(unique.celltype.atac) == 2){
    tmp.expand <- expand.grid(combn.celltypes[[1]], combn.celltypes[[2]])
    return(apply(tmp.expand[,1:2], 1, paste, collapse = "##"))
  }
  
  
}
