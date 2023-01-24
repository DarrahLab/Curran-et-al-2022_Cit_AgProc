###Cit-Abundance Relationships
#We are provided citrullinated positions in proteins, along with protmap abundance data
#for said proteins. We want to understand whether citrullination influences 
#the abundance of peptides in protmap. 
#Calculating differential expression of residues along linear sequence.
#AAG 6 April 2022 
rm(list = ls())
SAVE_FIGS <-FALSE
PARENT_DIR <- c("cit_abd") #This will be a subdir of master/results/
savedir <- paste0(getwd(), "/results/", PARENT_DIR)
if (!dir.exists(savedir) && SAVE_FIGS) {dir.create(savedir, showWarnings = FALSE)}

###WORKSPACE, LIBRARIES#####
library(strex)
library(seqinr)
library(ggplot2)
library(mutoss)
#setwd("/media/aag7319/WDBlue/ZZZ_protmap")
#setwd("E:/ZZZ_protmap")
#setwd("C:/Users/aag73/Desktop/ZZZ_protmap")
source('./functions/cit_corr_fn.r')
###LOAD DATA####
load("./objects/cit-sites-all-10.rda")
load("./objects/protmap_dat.rda")
###CALCULATE DIFFERENTIAL EXPRESSION OF RESIDUES IN MS DATA####
cutoff <- log(2, base = 2) # 100% or greater DiffEx
#Generate DiffEx matrices
Ags <- names(protmap_dat)
Ag_len <- c(353,866,491, 453, 466, 663)
diffEx <- vector(mode = "list", length = 6);  names(diffEx) <- Ags
for( a in 1:length(Ags)){
  dat <- protmap_dat[[a]]
  tmp <- matrix(data= NA, nrow = Ag_len[a], ncol = 3)
  colnames(tmp) <- c("pad2.native", "pad4.native", "pad4.pad2")
  
  n <- dat[[1]]; p2 <- dat[[2]]; 
  

  p2.n <- process1(p2,n)
 
  
  tmp[,1] <- p2.n

  if(a < 6){ #Don't run these is Ag is pad4
    p4 <- dat[[3]]
    p4.n <- process1(p4,n)
    p4.p2 <- process1(p4, p2)
    tmp[,2] <- p4.n
    tmp[,3] <- p4.p2
  }
  
  diffEx[[a]] <- tmp
  
}
rm(dat,tmp,n,p2,p4,p2.n,p4.n,p4.p2)
#Isolate residues of interest 
Up_Residues <- vector(mode = "list", length = 6);  names(Up_Residues) <- Ags
for( a in 1:length(Ags)){
  dat <- protmap_dat[[a]]
  tmp <- list(pad2.native = NA, pad4.native = NA, pad4.pad2 = NA) 
  tmp[[1]] <- which(diffEx[[a]][,1] >= cutoff)
  

  
  if(a < 6){ #Don't run these is Ag is pad4
    tmp[[2]]<- which(diffEx[[a]][,2] >= cutoff)
    tmp[[3]] <- which(diffEx[[a]][,3] >= cutoff)
  }
  
  Up_Residues[[a]] <- tmp
  
}
Down_Residues <- vector(mode = "list", length = 6);  names(Down_Residues) <- Ags
for( a in 1:length(Ags)){
  dat <- protmap_dat[[a]]
  tmp <- list(pad2.native = NA, pad4.native = NA, pad4.pad2 = NA) 
  tmp[[1]] <- which(diffEx[[a]][,1] < -cutoff)
  
  if(a < 6){ #Don't run these is Ag is pad4
    tmp[[2]]<- which(diffEx[[a]][,2] < -cutoff)
    tmp[[3]] <-which(diffEx[[a]][,2] < -cutoff)
  }
  
  Down_Residues[[a]] <- tmp
  
}
#save(diffEx, file = "./objects/diffExdat.rda")
###VISUALIZE DIFFEX####
p <- plot1(diffEx[[1]][,1], cutoff )

#Plot regular ol' Abundances 
if(SAVE_FIGS){
  
  subdir <- paste0(savedir, "/Abundances/")
  if (!dir.exists(subdir)){dir.create(subdir)}
  
  conds <- c("native", "pad2", "pad4")
  cutoff <- 6 #1e6
  for(f in 1: length(Ags)){
  c <- conds[1] #native
  a <- Ags[f]
  filename <- paste0(a,"_",c,"_abd.png")
  fullname <- paste0(subdir,filename)
  png(file = fullname, width = 750, height = 270)
  
  p <- plot0(protmap_dat[[f]][[1]], cutoff, paste(a,c))
  print(p)
  dev.off()
  
  
  c <- conds[2] #pad2
  a <- Ags[f]
  filename <- paste0(a,"_",c,"_abd.png")
  fullname <- paste0(subdir,filename)
  png(file = fullname, width = 750, height = 270)
  
  p <- plot0(protmap_dat[[f]][[2]],cutoff, paste(a,c))
  print(p)
  dev.off()
  
  if(f < 6){ #if not plotting pad4
  c <- conds[3] #pad4
  a <- Ags[f]
  filename <- paste0(a,"_",c,"_abd.png")
  fullname <- paste0(subdir,filename)
  png(file = fullname, width = 750, height = 270)
  
  p <- plot0(protmap_dat[[f]][[3]],cutoff, paste(a,c))
  print(p)
  dev.off()
  }
  }
  
  
}
#Plot Abd with Cit sites marked
if(SAVE_FIGS){
  
  subdir <- paste0(savedir, "/Abd_cit/")
  if (!dir.exists(subdir)){dir.create(subdir)}
  
  conds <- c("pad2", "pad4")
  cutoff <- 6 #1e6 unadj abundance
  
  Ags <- names(cit_sites_all[[1]])
  for(f in 1: length(Ags)){
    exg <- c(1:3,5,6)
    
    c <- conds[1] #p2.native
    a <- Ags[f]
    filename <- paste0(a,"_",c,"_abd.png")
    fullname <- paste0(subdir,filename)
    png(file = fullname, width = 750, height = 270)
    
    p <- plot3(protmap_dat[[exg[f]]][[2]], cit_sites_all[[1]][[f]], NULL, paste(a,c))
    print(p)
    dev.off()
    
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
  
  
}
#Plot DiffEx
if(SAVE_FIGS){
  
  subdir <- paste0(savedir, "/DiffEx/")
  if (!dir.exists(subdir)){dir.create(subdir)}
  
  conds <- c("pad2_native", "pad4_native", "pad4_pad2")
  cutoff <- log(2, base = 2) #50% diffEx
  Ags <- names(protmap_dat)
  for(f in 1: length(Ags)){
    c <- conds[1] #native
    a <- Ags[f]
    filename <- paste0(a,"_",c,"_abd.png")
    fullname <- paste0(subdir,filename)
    png(file = fullname, width = 750, height = 270)
    
    p <- plot1(diffEx[[f]][,1], cutoff, paste(a,c))
    print(p)
    dev.off()
    
    if(f < 6){
    c <- conds[2] #pad2
    a <- Ags[f]
    filename <- paste0(a,"_",c,"_abd.png")
    fullname <- paste0(subdir,filename)
    png(file = fullname, width = 750, height = 270)
    
    p <- plot1(diffEx[[f]][,2],cutoff, paste(a,c))
    print(p)
    dev.off()
    
     #if not plotting pad4
      c <- conds[3] #pad4
      a <- Ags[f]
      filename <- paste0(a,"_",c,"_abd.png")
      fullname <- paste0(subdir,filename)
      png(file = fullname, width = 750, height = 270)
      
      p <- plot1(diffEx[[f]][,3],cutoff, paste(a,c))
      print(p)
      dev.off()
    }
  }
  
  
}
#Plot DiffEx with Cit sites marked
if(SAVE_FIGS){
  
  subdir <- paste0(savedir, "/DiffEx_cit/")
  if (!dir.exists(subdir)){dir.create(subdir)}
  
  conds <- c("pad2_native", "pad4_native")
  cutoff <- log(2, base = 2) #100% diffEx
  
  Ags <- names(cit_sites_all[[1]])
  for(f in 1: length(Ags)){
    c <- conds[1] #p2.native
    a <- Ags[f]
    i <- c(1,2,3,5,6) #skip fib.g
    filename <- paste0(a,"_",c,"_abd.png")
    fullname <- paste0(subdir,filename)
    png(file = fullname, width = 750, height = 270)
    
    p <- plot41(diffEx[[i[f]]][,1], cit_sites_all[[1]][[f]], paste(a,c))
    print(p)
    dev.off()
    
    if(f < 5){
      c <- conds[2] #p4.native. Don't do if Pad4 is Ag
      a <- Ags[f]
      filename <- paste0(a,"_",c,"_abd.png")
      fullname <- paste0(subdir,filename)
      png(file = fullname, width = 750, height = 270)
      
      p <- plot41(diffEx[[f]][,2], cit_sites_all[[1]][[f]], cutoff, paste(a,c))
      print(p)
      dev.off()
      
  
    }
  }
  
  
}
