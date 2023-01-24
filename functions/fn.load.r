fn.load <- function(){
#Load all files in the './functions/' subdirectory without individually calling each.
#This is a function built for the lazy-ass. 
fn.list <- list.files('./functions/')
paths <- fn.list #placeholder
  for (i in 1: length(fn.list)){
    if (fn.list[i] != "fn.load.r"){
    paths[i] <- paste('./functions/', fn.list[i], sep = "")
    source(c(paths[i]))
    print(paste('Successfully added function ', fn.list[i], ' to environment.', sep= ""))
    }
  }
}