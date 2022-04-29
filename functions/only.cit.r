only.cit <- function(df.in){
  #AAG 7 November 2021
  #Expunge non-citrullination modifications from the modifications log
  #Doing this so that e.g. phospho peptide and peptide will be marked as 
  #identical and combined in downstream analysis. 
  #Not a perfect fn. If multiple mods are present and one of them is 
  #citrullination, non-cit changes will not be removed from record. 
  
  mod.idx <-  which((colnames(df.in) == "Modifications") == TRUE) 
  mods <-df.in[,mod.idx]
  deamid <- grep("Deamidated", mods)
  none <- which(is.na(mods))
  
  rm.idx <- c(deamid,none)
  
  df.out <- df.in
  df.out$Modifications[-rm.idx] <- NA
  
  return(df.out)
}