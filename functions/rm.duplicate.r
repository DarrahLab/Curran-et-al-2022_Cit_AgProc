rm.duplicate <- function(df.in, abd.idx){
  #AAG 7 Nov 2021
  #Combine duplicate peptides and delete extra row within dataframe. 
  #Determined duplicate if sequence and modifications are identical. 
  #Also leaving UID metadata column in output df. 
  
  df.out <- df.in
  mod.idx <-  which((colnames(df.in) == "Modifications") == TRUE) 
  as.idx <-  which((colnames(df.in) == "Annotated.Sequence") == TRUE) 
  
  mods <- df.in[,mod.idx]
  as <- df.in[,as.idx]
  df.out$UID <- paste(as,mods, sep = "::")
  nc <- ncol(df.out)
  df.out <- df.out[,c(nc, 1:(nc-1))]
  duplicates <- names(which(table(df.out$UID) > 1))
  n.dup <- length(duplicates)
  
  for (i in 1:n.dup){
    pos <- which(df.out$UID == duplicates[i])
    n <- length(pos)
    dat <- df.in[pos, abd.idx]
    dat <- as.matrix(sapply(dat, as.numeric))
    df.out[pos[1],(abd.idx+1)] <- colSums(dat)
    df.out <- df.out[-pos[2:n],]
  }

  return(df.out)
}