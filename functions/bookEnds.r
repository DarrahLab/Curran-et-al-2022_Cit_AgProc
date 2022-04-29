bookEnds <- function(pos.in, acc){
  #Accepts char vector of protmap 'Positions in Master Proteins'
  #returns matrix containing separated start and ends. 
  
  pos.spl <- unlist(strsplit(pos.in, ";"))
  if(!missing(acc)){pos.pr <- grep(acc, pos.spl)}else{pos.pr <- c(1:length(pos.spl))}
  st.sp <- read.table(text = pos.spl[pos.pr], fill = TRUE)[[2]]
  starts <- str_nth_number(st.sp, 1)
  stops <- str_nth_number(st.sp, 2)
  out  <- matrix(data = NA, nrow = length(pos.in), ncol = 2)
  out[,1] <- starts; out[,2] <- stops
  colnames(out) <- c("Start", "End")
  return(out)
}