concatenate.abundance <- function(master, addition, abd.idx){
  #AAG 7 November 2021
  #Combine df which may have unique peptides. 
  #For peptides in both df, simply add additional columns for new observations
  #For peptides only in 'addition', append new rows on 'master' to accomodate. 
  #Unique peptides in one df will have abundance NA in the other df.
  #Requires prior UID assignment and duplicate removal. 
  #UIDs must be in column 1 of both dataframes. 

n.rep <- length(abd.idx)
mas.in.add <- master[,1] %in% addition[,1]#Check UID equivalence 
add.in.mas <- addition[,1] %in% master[,1]
both.idx <- which(mas.in.add == TRUE)
mas.only.idx <- which(mas.in.add == FALSE)
add.only.idx <- which(add.in.mas == FALSE)
n.add <- length(add.only.idx)
n.mas <- length(mas.only.idx)
n.both <- length(both.idx)
new.names <- colnames(addition)[abd.idx]

#Check whether any unique peptides are in either sample. 
if (n.add == 0){ SKIP <- 1}
if(n.add != 0){SKIP <- 0}

#Increase size of master df to accomodate addition. 
if (!SKIP){
    #increase rows by #unique new peptides to add
    placehold <- data.frame(matrix(nrow = n.add, ncol = ncol(master)))
    colnames(placehold) <- colnames(master)
    master <- rbind(master,placehold)
}

#increase cols by #observations to add
placehold <- data.frame(matrix(ncol = n.rep, nrow = nrow(master)))
colnames(placehold) <- new.names
master <- cbind(master, placehold)

#Add peptides present in both groups
to.add <- which((add.in.mas == TRUE)) 
master[both.idx, (ncol(master) - n.rep + 1):ncol(master)] <- addition[to.add, abd.idx]

if (!SKIP){
#Add peptides only in addition to master
to.add <- which(add.in.mas == FALSE)
master[(nrow(master) - n.add + 1): nrow(master), -abd.idx] <- addition[to.add,]
}
return(master)
}