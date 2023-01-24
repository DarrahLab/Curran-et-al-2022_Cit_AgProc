##Various functions related to plotting pad protmap abundance
##Too lazy to cut/paste ggplot script all over the place
## AAG 9 November 2021 

#Plot single abd vector
plot.abd <- function(y, title,cutoff){
  #Edited 20 Feb 22
  if(missing(title)){title = NULL}
  if(missing(cutoff)){cutoff = NULL}

  id <- assign.id(list(y))
   
  
  if(!is.null(cutoff)){
   p <- ggplot(mapping = aes(x,y)) + geom_area(data = data.frame(x = c(1:length(y)), y = (y)))+ 
     labs(x =      paste(id, "AA Position", sep = " "), y = "Log2(Abundance)", title = title) + 
     theme_classic() + theme( axis.ticks.y = element_blank(), axis.text.y =element_blank(),
     panel.grid.major = element_line(colour = "light grey")) + 
     geom_hline(yintercept = cutoff_adj, color = "red")
   return(p)
  }else{
   p <-  ggplot(mapping = aes(x,y)) + geom_area(data = data.frame(x = c(1:length(y)), y = (y)))+ 
      labs(x =      paste(id, "AA Position", sep = " "), y = "Log2(Abundance)", title = title) + 
      theme_classic() + theme( axis.ticks.y = element_blank(), axis.text.y =element_blank(),
                               panel.grid.major = element_line(colour = "light grey"))
   return(p)
  }
}
plot.abd.inv <- function(y1,y2, title,cutoff){
  #Edited 20 Feb 22
  if(missing(title)){title = NULL}
  if(missing(cutoff)){cutoff = 2}
  if(length(y1) != length(y2)){
    print("Error: Provided vectors are unequal length") 
    return(NULL)
  }
  
  y1.zeros <- which(y1 == 0)
  y2.zeros <- which(y2 == 0)
  y1.unopposed <- y1.zeros[which((y1.zeros %in% y2.zeros) == FALSE)] #zero in y1 not in y2
  
  y2.unopposed <- y2.zeros[which((y2.zeros %in% y1.zeros) == FALSE)]
  common.zeros <- y1.zeros[which((y1.zeros %in% y2.zeros) == TRUE)]
  infs <- y1.unopposed
  neg.infs <- y2.unopposed 
  
  y <- y2 / y1
  biggest <-  max(abs(y))
  id <- assign.id(list(y))
  y[common.zeros] = 0
  y[y1.unopposed] = biggest + 3 #  + Inf 
  y[y2.unopposed] = -biggest - 3 # - Inf

  
  
  
  ggplot(mapping = aes(x,y)) + geom_area(data = data.frame(x = c(1:length(y)), y = (y)))+ 
    labs(x =      paste(id, "AA Position", sep = " "), y = "Log2(Abundance)", title = title) + 
    theme_classic() + theme( axis.ticks.y = element_blank(), axis.text.y =element_blank(),
    panel.grid.major = element_line(colour = "light grey")) + 
    geom_hline(yintercept = cutoff, color = "red") + 
    geom_hline(yintercept = -cutoff, color = "red")
}

#When given nx3 matrix with native, pad2, pad4 abd, plot native and pad2
#if inv = 1 do relative enrichment (subtract the inputs)
plot.abd.pad2 <- function(y, inv,lim){
  id <- assign.id(y)
if (!inv){
  y[,1] <- log2t(y[,1])
  y[,2] <- log2t(y[,2])
  return(ggplot(mapping = aes(x,y)) + 
    geom_area(data = data.frame(x = c(1:dim(y)[1]), y = (y[,1])) , width = 0.4, stat = 'identity', alpha = 0.4, fill = "black") +
    geom_area(data = data.frame(x = c(1:dim(y)[1]), y = (y[,2])) , width = 0.4, stat = 'identity', alpha = 0.4, fill = "dark blue") +
    theme_classic() + scale_y_continuous(expand = c(0,0)))
    #labs(x = paste(id, "AA Position", sep = " "),  y = "Log2 Abundance", title = c("Native and PAD2")))
  
}
if(inv){
  y <- log2t(y)
  y = ((y[,2] - y[,1]))
  y[y>lim[2]] <- lim[2]
  y[y<lim[1]] <- lim[1]
  return(ggplot(mapping = aes(x,y)) + 
    geom_area(data = data.frame(x = c(1:length(y)), y = y) , width = 0.4, stat = 'identity', alpha = 1, fill = "black") +
    scale_y_continuous(expand = c(0,0), limits = lim) + geom_hline(yintercept = 0) + 
      theme_classic() +  theme(text = element_text(size = 15)) + 
     labs(x = "",  y = "", title = ""))
}
  
}

#When given nx3 matrix with native, pad2, pad4 abd, plot native and pad4
#if inv = 1 do relative enrichment (subtract the inputs)
plot.abd.pad4 <- function(y, inv,lim){
  id <- assign.id(y)
  if (!inv){
    y[,1] <- log2t(y[,1])
    y[,3] <- log2t(y[,3])
    return(ggplot(mapping = aes(x,y)) + 
             geom_area(data = data.frame(x = c(1:dim(y)[1]), y = (y[,1])) , width = 0.4, stat = 'identity', alpha = 0.4, fill = "black") +
             geom_area(data = data.frame(x = c(1:dim(y)[1]), y = (y[,3])) , width = 0.4, stat = 'identity', alpha = 0.4, fill = "dark red") +
             theme_classic() + scale_y_continuous(expand = c(0,0)))  
            # labs(x = paste(id, "AA Position", sep = " "),  y = "Log2 Abundance", title = c("Native and PAD4")))
    
  }
  if(inv){
    y <- log2t(y)
    y = ((y[,3] - y[,1]))
    y[y>lim[2]] <- lim[2]
    y[y<lim[1]] <- lim[1]
    return(ggplot(mapping = aes(x,y)) + 
             geom_area(data = data.frame(x = c(1:length(y)), y = y) , width = 0.4, stat = 'identity', alpha = 1, fill = "black") +
             scale_y_continuous(expand = c(0,0), limits = lim) + geom_hline(yintercept = 0) + 
             theme_classic() + theme(text = element_text(size = 15)) + 
             labs(x = "",  y = "", title = ""))
    
  }
}

#When given nx3 matrix with native, pad2, pad4 abd, plot pad2 and pad4
#if inv = 1 do relative enrichment (subtract the inputs)
plot.abd.pads <- function(y, inv){
  id <- assign.id(y)
  if (!inv){
    y[,2] <- log2t(y[,2])
    y[,3] <- log2t(y[,3])
    return(ggplot(mapping = aes(x,y)) + 
             geom_area(data = data.frame(x = c(1:dim(y)[1]), y = (y[,2])) , width = 0.4, stat = 'identity', alpha = 0.4, fill = "dark blue") +
             geom_area(data = data.frame(x = c(1:dim(y)[1]), y = (y[,3])) , width = 0.4, stat = 'identity', alpha = 0.4, fill = "dark red") +
             theme_classic() + scale_y_continuous(expand = c(0,0)) + 
             labs(x = paste(id, "AA Position", sep = " "),  y = "Log2 Abundance", title = c("PAD2 and PAD4")))
    
  }
  if(inv){
    y[,2] <- log2t(y[,2])
    y[,3] <- log2t(y[,3])
    return(ggplot(mapping = aes(x,y)) + 
             geom_area(data = data.frame(x = c(1:dim(y)[1]), y = ((y[,3] - y[,2]))) , width = 0.4, stat = 'identity', alpha = 1, fill = "black") +
             theme_classic() + scale_y_continuous(expand = c(0,0)) + geom_hline(yintercept = 0) + 
             labs(x = paste(id, "AA Position", sep = " "),  y = " Log2 Enrichment", title = c("PAD4 v. PAD2")))
    
  }
}

#If you really want to, plot all three (native, pad2, pad4) at once. 
#Do not especially recommend.
plot.abd.3 <- function(y, title){
  id <- assign.id(y)
  
  ggplot(mapping = aes(x,y)) + 
    geom_area(data = data.frame(x = c(1:dim(y)[1]), y = (y[,1])) , width = 0.4, stat = 'identity', alpha = 0.4, fill = "black") +
    geom_area(data = data.frame(x = c(1:dim(y)[1]), y = (y[,2])) , width = 0.4, stat = 'identity', alpha = 0.4, fill = "dark blue") +
    geom_area(data = data.frame(x = c(1:dim(y)[1]), y = (y[,3])) , width = 0.4, stat = 'identity', alpha = 0.4, fill = "dark green") +
    theme_classic() + scale_y_continuous(expand = c(0,0), trans = 'log2') + 
    labs(x = paste(id, "AA Position", sep = " "),  y = "Log2 Abundance", title = title)
}

#Find ID of input abundance matrix by AA length. 
assign.id <- function(y){
  if (length(y[[1]]) == 888){ id <- c("HMGCR")} else
    if (length(y[[1]]) == 866){ id <- c("Fib. Alpha")} else 
      if (length(y[[1]]) == 491){ id <- c("Fib. Beta")} else 
        if (length(y[[1]]) == 453){ id <- c("Fib. Gamma")}else
          if (length(y[[1]]) == 353){ id <- c("RA33")} else 
            if (length(y[[1]]) == 466){ id <- c("Vim")} else 
              if (length(y[[1]]) == 663){ id <- c("Pad4")} else {
                print("Error: Input length not recognized as known protein")
                break}
    return(id)
}
#Log2 transform data prior to plotting. 
log2t<- function(y){
  out <- y + 1
  out <- log(out, base = 2)
  return(out)
}