#Final protmap cityscape plots. 
#AAG 6 Apr 2022
setwd("C:/Users/aag73/Desktop/ZZZ_protmap")

plot42 <- function(d, cit_sites, cutoff, specifier){
  if(missing(cutoff)){
    cutoff <- NULL
  }
  if(missing(specifier)){
    specifier <- NULL
  }
  t <- paste0("Diff Ex: ", specifier)
  #lim <- max(abs(d[abs(d) < Inf]))
  lim <- 9.5
  
  #Make cit_sites df 
  x <- cit_sites
  y <- 0
  cits <- data.frame(x = x, y = y)
  cits2 <- data.frame(x = x, y = -y)
  
  
  p <- ggplot(data = data.frame(x = c(1:length(d)), y = d), aes(x,y)) + 
    geom_area( color = "black", fill = "black") + 
    labs(title = t, y = "log2(Abd)") + theme_minimal() + theme(text = element_text(size = 15)) + 
    geom_hline(yintercept = 0, color = "black")  + 
    #theme_prism() + 
    geom_point(data = cits, col = "black", fill = "#00AD66", shape =21,  size = 4) + 
    scale_x_continuous(breaks = seq(1,350, by = 20) - 1, minor_breaks = seq(1,350, by = 5) - 1) + 
    #ylim(-lim,lim)
    ylim(-5,5)
  
  return(p)
  
  
} 
###WORKSPACE, LIBRARIES#####
library(strex)
library(seqinr)
library(ggplot2)
library(mutoss)
#setwd("/media/aag7319/WDBlue/ZZZ_protmap")
#setwd("E:/ZZZ_protmap")
setwd("C:/Users/aag73/Desktop/ZZZ_protmap")

###LOAD DATA####
load("./objects/cit-sites-all-10.rda")
load("./objects/protmap_dat.rda")
load("./objects/diffExdat.rda")

###VISUALIZE DIFFEX####
#Plot DiffEx with Cit sites marked
  try = "grBl2"
  subdir <- paste0(savedir, "/DiffEx_cit/")
  if (!dir.exists(subdir)){dir.create(subdir)}
  
  conds <- c("pad2_native", "pad4_native")
  cutoff <- log(2, base = 2) #100% diffEx
  
  Ags <- names(cit_sites_all[[1]])
  for(f in 1:length(Ags)){
    c <- conds[1] #p2.native
    a <- Ags[f]
    i <- c(1,2,3,5,6) #skip fib.g
    filename <- paste0(a,"_",c,"_",try, "_abd.png")
    fullname <- paste0(subdir,filename)
    png(file = fullname, width = 750, height = 270)
    
    p <- plot42(diffEx[[i[f]]][,1], cit_sites_all[[1]][[f]], paste(a,c))
    print(p)
   
    
    if(f < 5){
      c <- conds[2] #p4.native. Don't do if Pad4 is Ag
      a <- Ags[f]
      filename <- paste0(a,"_",c,"_abd.png")
      fullname <- paste0(subdir,filename)
     
      p <- plot42(diffEx[[f]][,2], cit_sites_all[[1]][[f]], cutoff, paste(a,c))
      print(p)
   
      
      
    }
  }
  
  


#Plot regular ol' Abundances 

  
  subdir <- paste0(savedir, "/Abundances/")
 
  
  conds <- c("native", "pad2", "pad4")
  cutoff <- 6 #1e6
  for(f in 1: length(Ags)){
    c <- conds[1] #native
    a <- Ags[f]
    filename <- paste0(a,"_",c,"_abd.png")
    fullname <- paste0(subdir,filename)

    
    p <- plot0(protmap_dat[[f]][[1]], cutoff, paste(a,c))
    print(p)

    
    
    c <- conds[2] #pad2
    a <- Ags[f]
    filename <- paste0(a,"_",c,"_abd.png")
    fullname <- paste0(subdir,filename)

    
    p <- plot0(protmap_dat[[f]][[2]],cutoff, paste(a,c))
    print(p)

    
    if(f < 6){ #if not plotting pad4
      c <- conds[3] #pad4
      a <- Ags[f]
      filename <- paste0(a,"_",c,"_abd.png")
      fullname <- paste0(subdir,filename)
 
      
      p <- plot0(protmap_dat[[f]][[3]],cutoff, paste(a,c))
      print(p)

    }
  }
  
  

#Plot Abd with Cit sites marked

  
  subdir <- paste0(savedir, "/Abd_cit/")
 
  
  conds <- c("pad2", "pad4")
  cutoff <- 6 #1e6 unadj abundance
  
  Ags <- names(cit_sites_all[[1]])
  for(f in 1: length(Ags)){
    exg <- c(1:3,5,6)
    
    c <- conds[1] #p2.native
    a <- Ags[f]
    filename <- paste0(a,"_",c,"_abd.png")
    fullname <- paste0(subdir,filename)
 
    
    p <- plot3(protmap_dat[[exg[f]]][[2]], cit_sites_all[[1]][[f]], NULL, paste(a,c))

    
    if(f < 5){
      c <- conds[2] #p4.native. Don't do if Pad4 is Ag
      a <- Ags[f]
      filename <- paste0(a,"_",c,"_abd.png")
      fullname <- paste0(subdir,filename)
      png(file = fullname, width = 750, height = 270)
      
      p <- plot3(protmap_dat[[exg[f]]][[2]], cit_sites_all[[1]][[f]], NULL, paste(a,c))
      print(p)
      dev.off()
      
      
    }
  }
  
  

#Plot DiffEx

  
  subdir <- paste0(savedir, "/DiffEx/")
 
  
  conds <- c("pad2_native", "pad4_native", "pad4_pad2")
  cutoff <- log(2, base = 2) #50% diffEx
  Ags <- names(protmap_dat)
  for(f in 1: length(Ags)){
    c <- conds[1] #native
    a <- Ags[f]
    filename <- paste0(a,"_",c,"_abd.png")
    fullname <- paste0(subdir,filename)
  
    p <- plot42(diffEx[[f]][,1], cutoff, paste(a,c))
    print(p)

    
    if(f < 6){
      c <- conds[2] #pad2
      a <- Ags[f]
      filename <- paste0(a,"_",c,"_abd.png")
      fullname <- paste0(subdir,filename)
      
      
      p <- plot42(diffEx[[f]][,2],cutoff, paste(a,c))
      print(p)

      
      #if not plotting pad4
      c <- conds[3] #pad4
      a <- Ags[f]
      filename <- paste0(a,"_",c,"_abd.png")
      fullname <- paste0(subdir,filename)
     
      p <- plot42(diffEx[[f]][,3],cutoff, paste(a,c))
      print(p)

    }
  }
  


