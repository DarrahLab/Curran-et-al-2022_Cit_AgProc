find_deamidation <- function(df.in){
  #I may have originally written this fn to take multiple prot fragments at once 
  #(ie dim(df.in)[1] > 1) but that is not how I am using it in larger script anymore. 
  #So I am modifying line 18 from df.in$start[l] to having one start var. 
 cit.idx <- grep("Deamidated", df.in$Modifications)
 cit <- df.in$Modifications[cit.idx]
 n.pep <- length(cit.idx)
 
 tmp <- sub(".*Deamidated ", "",cit)
 tmp <- strsplit(tmp, split = ";")
 
#Identify non-deamidation modifications on this same fragment, and filter them. 
 to.rm <- c( grep("Oxidation", tmp[[1]]), grep("Carbamidomethyl", tmp[[1]]))
          
 if(length(to.rm) > 0) {tmp <- tmp[[1]][-to.rm]}else{tmp <- tmp[[1]]}
 
 to.keep <- grep("R", tmp)
 tmp <- tmp[to.keep]
 tmp <- lapply(tmp, str_first_number)

 for (l in 1:length(tmp)){
   tmp[[l]] <- tmp[[l]] + df.in$Start - 1
   #tmp[[l]] <- tmp[[l]] + df.in$Start[l] - 1
 }
 cit.global.sites <- unlist(tmp)
 #rm citrullinations w unspecified location
 cit.global.sites <- cit.global.sites[which(!is.na(cit.global.sites))]
 
 return(unique(cit.global.sites))
}