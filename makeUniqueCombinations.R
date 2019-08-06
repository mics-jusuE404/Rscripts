## Given two lists of elements (with unique names), print all pairwise and fully unique combinations:
makeUniqueCombinations <- function(list.a, list.b){
  
  ## in case of trouble, stop without this ugly stop() error message,
  ## credit: https://stackoverflow.com/questions/14469522/stop-an-r-program-without-error
  stopQuietly <- function() {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }
  
  if (0 %in% sapply(c(list.a, list.b), function(voltorb) nchar(voltorb))) {
    message("[Error] Check the input lists, seems that at least one contains empty strings or other odd entries")
    stopQuietly()
  }
  
  ## put some unique letters around the input names to ensure they are unique:
  list.a <- paste0("makeAunique", list.a, "makeAunique")
  list.b <- paste0("makeBunique", list.b, "makeBunique")
  
  ## make sure elements are unique in each list:
  if ( length(list.a) != length(unique(list.a)) | 
       length(list.b) != length(unique(list.b)) ){
    message("[Error]: At least one of the two lists contains duplicated entries!")
    stopQuietly()
  }
  
  ## make all possible pairwise combinations:
  combos.all <- outer(list.a,list.b, FUN=paste)
  
  ## separate by "__" as this is pretty sure a unique delimiter
  combos.all <- as.character(sapply(c(combos.all), function(x) gsub(" ", "__",x)))
  
  ## if less than three combinations (so only 1 element in either list.a or list.b),
  ## return data at this point as combos.all is already the full list of valid combinations:
  if (length(combos.all) < 3){
    message(paste("[Note]: At least one of both lists contains only a single entry.",
                  paste("        Returning", length(combos.all), "unique combinations"), 
                  sep="\n")
    )
    message("        The first part of the string should be the ATAC-seq sample, the second one the RNA-seq.")
    combos.all <- gsub("makeAunique|makeBunique", "", paste0("makeAunique", combos.all))
    return(combos.all)
  } 
  
  ## Now for each of combos.all find unqieu combinations with all other elements,
  ## without that any of the single elements is duplicated:
  combos.unique <- sapply(combos.all, function(x){
    
    ## first remove from combos.all all elements that overlap
    tmp.combn <- combos.all[!grepl(gsub("__","|", x), combos.all)]
    
    ## now make sure that partners are compatible among each other:
    return(
      sapply(tmp.combn, function(k) {
      return(
        sort(paste(x, k, tmp.combn[!grepl(gsub("__","|",k), tmp.combn)]))
        )
      })
    )
  })
  
  tmp.out <- sapply(combos.unique, function(f){
    paste(sort(strsplit(f, split=" ")[[1]]), collapse = "@@")
    }
  )
  
  tmp.out <- unique(tmp.out)
  
  tmp.out <- gsub("makeAunique|makeBunique", "", paste0("makeAunique", tmp.out))
  
  message(paste("[Info]: Returning", length(tmp.out), "unique combinations."))
  message("        The first part of the string should be the ATAC-seq sample, the second one the RNA-seq.")
  
  return(tmp.out)
  
}

## the names (which will be the colnames of the rna and atac data will then be used to make the inputs for InTAD)
list.a <- c("sample1_ATAC", "sampe2_ATAC", "sample3_ATAC")
list.b <- c("sample1_RNA", "sample2_RNA", "sample3_RNA")

makeUniqueCombinations(list.a = list.a, list.b = list.b)
