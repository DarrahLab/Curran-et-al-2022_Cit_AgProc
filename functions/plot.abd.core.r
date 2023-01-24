##Various functions related to plotting pad protmap abundance
##Too lazy to cut/paste ggplot script all over the place
## AAG 9 November 2021 

#Plot single abd vector
plot.abd <- function(y, title){
  if (length(y) == 888){ id <- c("HMGCR")} else
    if (length(y) == 866){ id <- c("Fib. Alpha")} else 
      if (length(y) == 491){ id <- c("Fib. Beta")} else 
        if (length(y) == 453){ id <- c("Fib. Gamma")}else {id <- ""}
  ggplot(mapping = aes(x,y)) + 
    geom_area(data = data.frame(x = c(1:length(y)), y = (y)) , width = 0.4, stat = 'identity') +
    theme_classic() + scale_y_continuous(expand = c(0,0), trans = 'log2') + labs(x = paste(id, "AA Position", sep = " "), 
                                                                                 y = "Abundance", title = title)
}

#When given nx3 matrix with native, pad2, pad4 abd, plot native and pad2
#if inv = 1 do relative enrichment (subtract the inputs)
plot.abd.pad2.C <- function(y, c){
  if (dim(y)[1] == 866){ id <- c("Fib. Alpha"); z = y * c[[1]]; c = c[[1]]} else 
    if (dim(y)[1] == 491){ id <- c("Fib. Beta"); z = y * c[[2]]} else 
      if (dim(y)[1] == 453){ id <- c("Fib. Gamma"); z = y * c[[3]]}else {id <- ""}
  
  y[,1] <- log2t(y[,1])
  y[,2] <- log2t(y[,2])
  z <- log2t(z)
  
    return(ggplot(mapping = aes(x,y)) + 
             geom_area(data = data.frame(x = c(1:dim(y)[1]), y = ((y[,2] - y[,1]))) , width = 0.4, stat = 'identity', alpha = 1, fill = "black") +
             theme_classic() + scale_y_continuous(expand = c(0,0)) + geom_hline(yintercept = 0) + 
             geom_area(data = data.frame(x = c(1:dim(y)[1]), y = (z[,2] - z[,1])) , width = 0.4, stat = 'identity', alpha =1 , fill = "red") +
             labs(x = paste(id, "AA Position", sep = " "),  y = "Log2 Enrichment", title = c("PAD2 v. Native")))
    
  
}

#When given nx3 matrix with native, pad2, pad4 abd, plot native and pad4
#if inv = 1 do relative enrichment (subtract the inputs)
plot.abd.pad4.C <- function(y, c){
  if (dim(y)[1] == 866){ id <- c("Fib. Alpha"); z = y * c[[1]]; c = c[[1]]} else 
    if (dim(y)[1] == 491){ id <- c("Fib. Beta"); z = y * c[[2]]; c = c[[2]]} else 
      if (dim(y)[1] == 453){ id <- c("Fib. Gamma"); z = y * c[[3]]; c = c[[3]]}else {id <- ""}
  y[,1] <- log2t(y[,1])
  y[,3] <- log2t(y[,3])
  z <- log2t(z)

    return(ggplot(mapping = aes(x,y)) + 
             geom_area(data = data.frame(x = c(1:dim(y)[1]), y = ((y[,3] - y[,1]))) , width = 0.4, stat = 'identity', alpha = 1, fill = "black") +
             theme_classic() + scale_y_continuous(expand = c(0,0)) + geom_hline(yintercept = 0) + 
             geom_area(data = data.frame(x = c(1:dim(y)[1]), y = (z[,3] - z[,1])) , width = 0.4, stat = 'identity', alpha =1 , fill = "red") +
             labs(x = paste(id, "AA Position", sep = " "),  y = " Log2 Enrichment", title = c("PAD4 v. Native")))
    
  
}

#When given nx3 matrix with native, pad2, pad4 abd, plot pad2 and pad4
#if inv = 1 do relative enrichment (subtract the inputs)
plot.abd.pads.C <- function(y,c){
  if (dim(y)[1] == 866){ id <- c("Fib. Alpha"); z = y * c[[1]]; c = c[[1]]} else 
    if (dim(y)[1] == 491){ id <- c("Fib. Beta"); z = y * c[[2]]; c = c[[2]]} else 
      if (dim(y)[1] == 453){ id <- c("Fib. Gamma"); z = y * c[[3]]; c = c[[3]]}else {id <- ""}

    y[,2] <- log2t(y[,2])
    y[,3] <- log2t(y[,3])
    z <- log2t(z)
    return(ggplot(mapping = aes(x,y)) + 
             geom_area(data = data.frame(x = c(1:dim(y)[1]), y = ((y[,3] - y[,2]))) , width = 0.4, stat = 'identity', alpha = 1, fill = "black") +
             geom_area(data = data.frame(x = c(1:dim(y)[1]), y = (z[,3] - z[,2])) , width = 0.4, stat = 'identity', alpha =1 , fill = "red") +
             theme_classic() + scale_y_continuous(expand = c(0,0)) + geom_hline(yintercept = 0) + 
             labs(x = paste(id, "AA Position", sep = " "),  y = " Log2 Enrichment", title = c("PAD4 v. PAD2")))
    
}

#If you really want to, plot all three (native, pad2, pad4) at once. 
#Do not especially recommend.
plot.abd.3 <- function(y, title){
  if (dim(y)[1] == 866){ id <- c("Fib. Alpha")} else 
    if (dim(y)[1] == 491){ id <- c("Fib. Beta")} else 
      if (dim(y)[1] == 453){ id <- c("Fib. Gamma")}else {id <- ""}
  
  ggplot(mapping = aes(x,y)) + 
    geom_area(data = data.frame(x = c(1:dim(y)[1]), y = (y[,1])) , width = 0.4, stat = 'identity', alpha = 0.4, fill = "black") +
    geom_area(data = data.frame(x = c(1:dim(y)[1]), y = (y[,2])) , width = 0.4, stat = 'identity', alpha = 0.4, fill = "dark blue") +
    geom_area(data = data.frame(x = c(1:dim(y)[1]), y = (y[,3])) , width = 0.4, stat = 'identity', alpha = 0.4, fill = "dark green") +
    theme_classic() + scale_y_continuous(expand = c(0,0), trans = 'log2') + 
    labs(x = paste(id, "AA Position", sep = " "),  y = "Log2 Abundance", title = title)
}

#Log2 transform data prior to plotting. 
log2t<- function(y){
  out <- y + 1
  out <- log(out, base = 2)
  return(out)
}