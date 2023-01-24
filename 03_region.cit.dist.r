#Distance to cit No. 2
#Isolate 'Created' and 'Destroyed' regions of MS data
#For each region, find out how far it is from a citrulline. 
#Three modes: bi-directional, fwd, rev, +/- intra-cit 
# AAG 11 March 2022

rm(list = ls())
SAVE_FIGS <-FALSE
PARENT_DIR <- c("cit_abd") #This will be a subdir of master/results/
savedir <- paste0(getwd(), "/results/", PARENT_DIR)
if (!dir.exists(savedir) && SAVE_FIGS) {dir.create(savedir, showWarnings = FALSE)}

##Libraries / Workdir / load Data----
setwd("/media/aag7319/WDBlue/ZZZ_protmap")

library(ggplot2)
library(ggprism)
load("./objects/diffExdat.rda")
load("./objects/cit-sites-all-10.rda")
load("./objects/protmap_dat.rda")
diffEx <- diffEx[-4] #remove fib.gamma, which has no cit sites. 
protmap_dat <- protmap_dat[-4]
#Functions----
regionMake <- function(diffEx_vec, thresh = 1){
  #Return dataframe containing the start, end, len of each highly expressed region
  #DiffEx_vec log2 normalized
  len <- length(diffEx_vec)
  logiC <- diffEx_vec >= thresh
  logiD <- diffEx_vec <= -thresh
  
  logis <- list(logiC, logiD)
  rdf <- data.frame(start = NA, end = NA, type = NA)
  
  for ( log in 1:2){ #Perform twice for created and destroyed regions. 
  logi <- logis[[log]]
  diffy <- logi[2:len] - logi[1:(len - 1)]
  starts <- which(diffy == 1) + 1
  ends <- which(diffy == -1)
  type <- c("Created", "Destroyed")[log]
  
  nSt <- length(starts)
  nEn <- length(ends)
  #Check for 0 - E1
  if (ends[1] < starts[1]){
    tmp <- data.frame(start = 1, end = ends[1], type = type)
    rdf <- rbind(rdf, tmp)
    ends <- ends[-1]
    nEn <- length(ends)
  }
  if (ends[nEn] < starts[nSt]){
    tmp <- data.frame(start = starts[nSt], end = len, type = type)
    rdf <- rbind(rdf, tmp)
    starts <- starts[-nSt]
    nSt <- length(starts)
  }
  
  #After trimming exceptions, the remaining starts and ends 
  #should be paired
  tmp <- data.frame(start = starts, end = ends, type = type)
  rdf <- rbind(rdf,tmp)
  
  }
  
  rdf <- rdf[-1,]
  rdf$len <- rdf$end - rdf$start + 1
  o <- order(rdf$start)
  rdf <- rdf[o,]
  return(rdf)
}
regionDist <- function(rdf, cits, dir = 0, intra = TRUE){
  #returns vector containing distance to cit from each region in rdf
  #dir = 1 searches only forward, -1 only backward, 0 bidirectional
  #intra = TRUE means any cit within a region gives cit dist = 0. if False, 
  #intra-region cits are ignored in computing distance. 
  
  nR <- dim(rdf)[1]
  dists <- vector(mode = "integer", length = nR)
  for (r in 1:nR){
    st <- rdf[r,1]
    en <- rdf[r,2]
    #Check whether cit site within region. Return 0 if appropriate. 
    if(any(cits >= st & cits <= en) & intra == TRUE){
      dists[r] <- 0
      next 
    }
    
    if(dir ==1){ #fwd direction
      Cp <- cits[which(cits > en)[1]]
      if(length(Cp) == 0){
        dists[r] <- NA
        next
      }
      dists[r] <- Cp - en
      next
    }else if(dir == -1){ #rev direction
      Cp <- cits[which(cits < st)[length(which(cits<st))]]
      if(length(Cp) == 0){
        dists[r] <- NA
        next
      }
      dists[r] <- st - Cp
      next
    }else if(dir == 0){
      if (any(cits > en)){
      Cpf <- cits[which(cits > en)[1]]
      }else{Cpf = NA}
      if(any(cits < st)){
      Cpr <- cits[which(cits < st)[length(which(cits<st))]]
      }else{Cpr = NA}
      
      df <- Cpf - en
      dr <- st - Cpr
      
      if(is.na(Cpf)){
        dists[r] <- dr
      }else if(is.na(Cpr)){
        dists[r] <- df
      }else{
        dists[r] <- min(df,dr)
      }
    }
  }
  
  return(dists)
}
rcomb <- function(Y,N,M){
  
  #Y max value in choosing range 
  #N choose how many? 
  #M how many combinations do you want 
  #example rcomb(10,2,5) Give me 5 combinations of 10 choose 2
  set.seed(42)
  x <- seq(1, Y)
  out <- list()
  while(length(out) < M) {
    out <- c(out,
             unique(replicate(M - length(out), sort(sample(x, N)), simplify = FALSE)))
  }
  
  return(out)
  
}
diffHist <- function(diffres_vec){
  p <- ggplot(data = data.frame(x = diffres_vec), aes(x = x)) + 
      geom_histogram(bins = 50) + theme_classic()
  return(p)
}

diffDens <- function(diffres_vec){
  p <- ggplot(data = data.frame(x = diffres_vec), aes(x = x)) + 
    geom_density() + theme_classic()
  return(p)
}
diffVln <- function(diffres_vec){
  p <- ggplot(data = data.frame(x = diffres_vec), aes(x = x)) + 
    geom_violin() + theme_classic()
  return(p)
}

diffDens24 <-function(d2,d4){
  #black = d2, blue = d4
  p <- ggplot(data = data.frame(x = d2), aes(x = x)) + 
    geom_density(fill = "black", alpha = 0.6) + theme_classic() + 
    geom_density(data = data.frame(x = d4), fill = "blue", color = "blue", alpha = 0.6)
}
#Ascertain protein coverage - 'deep' reads####
#Find sites in pad2 and pad4 conditions abd > cutoff (1e6) We say these 
#are read at enough depth to call a citrulline if it twere citrullinated. 
Ags <- names(protmap_dat)
deep_reads <- vector(mode = "list", length = 5); names(deep_reads) <- Ags
cutoff <- 1e6
for(a in 1:length(Ags)){
  tmp <- list(pad2 = NA, pad4 = NA)
  
  p2.deep <- which(protmap_dat[[a]][[2]] >= cutoff)
  tmp[[1]] <- p2.deep
  
  if(a < 5){
    p4.deep <- which(protmap_dat[[a]][[3]] >= cutoff)
    tmp[[2]] <- p4.deep
  
  }
  
  deep_reads[[a]] <- tmp
}
rm(tmp, p2.deep, p4.deep)
#Distance to cit----
#Test parameters 
N <- 2000 #number of simulated cit sites 
thr <- 1 #threshold for 'created' region. Default = 1 (2x expressed)
intra <- TRUE #if cit inside region, is d-to-cit 0? If false, ignore inner cits. 

#placeholder for data. 
tmp <- vector(mode = "list", length = 5); names(tmp) <- names(diffEx)
diffres <- list(tmp,tmp); names(diffres) <- c("pad2", "pad4"); rm(tmp)
diffres <- list(diffres,diffres,diffres); names(diffres) <- c("bidir", "fwd", "rev")
pvals <- diffres
distances <- diffres
expectations <- diffres
brdf <- as.data.frame.matrix(matrix(data = NA, ncol = 9, nrow = 1))
names(brdf) <- c("start", "end", "type", "len", "ag", "cond", "dist", "exp", "diff")
#execute test
for (Ag in 1:5){ #for each antigen
  for(pad in 1:2){ #testing pad2 or pad4 cit sites 
    if(Ag == 5 & pad == 2){break}
    cits <- cit_sites_all[[pad]][[Ag]]
    diffVec <- diffEx[[Ag]][,pad]
    #seqlen <- length(diffVec)
    seqlen <- length(deep_reads[[Ag]][[pad]])
    ncits <- length(cits)
    bootstrap_cits <- rcomb(seqlen, ncits, N )
    
    rdf <- regionMake(diffVec, thresh = thr)
    for(d in 1){ #For each direction mode (actually just bidir. Set 1:3 for all settings. )
      dir = c(0,1,-1)[d]
      dists <- regionDist(rdf,cits,dir, intra)
      
      d_exp <- matrix(data = NA, nrow = N, ncol = dim(rdf)[1])
      for(b in 1:N){ #Run bootstrap simulation
        seqlen <- length(diffVec)
        boocits <- deep_reads[[Ag]][[pad]][bootstrap_cits[[b]]]
        d_exp[b,] <- regionDist(rdf, bootstrap_cits[[b]], dir)
      }
   
      #Get mean expected values. 
      exp <- colSums(d_exp, na.rm = TRUE) / N  #Think about NAs. Unsure what to do. 
      diffres[[d]][[pad]][[Ag]] <- dists - exp
      distances[[d]][[pad]][[Ag]] <- dists 
      expectations[[d]][[pad]][[Ag]] <-  exp
      rdf$ag <- names(diffEx)[Ag]
      rdf$cond <- c("pad2", "pad4")[pad] 
      rdf$dist <- dists
      rdf$exp <- exp
      rdf$diff <- dists-exp
      
      brdf <- rbind(brdf,rdf)
      #Get p values
      ps <- vector(mode = "double", length = dim(rdf)[1])
      for(r in 1:length(dists)){
      #Find number of d_exp less than actual.
        p <- sum(d_exp[,r] <= dists[r], na.rm = TRUE) / N
        ps[r] <- min(p, 1-p)
      }
      pvals[[d]][[pad]][[Ag]] <- ps
    }
  }
}
brdf <- brdf[-1,]
p.idx <- which(brdf$ag == "pad4")
brdf$cond[p.idx] <- "pad2"

#Fig A. Bargraph of Created and Destroyed regions w and wo cits. ----
to.keep <- match(c("type", "ag", "dist", "cond"), colnames(brdf))
df1 <- brdf[,to.keep]
df1$cit <- as.factor(!as.logical(df1$dist))#if cit within region, TRUE
df1$dist <- NULL
df1$cond <- as.factor(df1$cond)

df1$cond[which(df1$ag == "pad4")] <- "pad4"


ggplot(df1, aes(fill = cit, x = type)) + 
  geom_bar(stat = "count") + 
  theme_prism() + 
  xlab(NULL) + 
  ylab("Count") + 
  labs(title = "Cits in DiffEx Regions") + 
  #ylim(0,60) + 
  facet_grid(~ cond)
  


ggplot(df1[-which(df1$ag == "pad4"),], aes(fill = cit, x = type)) + 
  geom_bar(stat = "count") + 
  theme_prism() + 
  xlab(NULL) + 
  ylab("Count") + 
  labs(title = "Cits in DiffEx Regions") + 
  
  facet_grid(~ cond)





  ggplot(df2, aes(fill = cit, x = type)) + 
  geom_bar(stat = "count") + 
  theme_prism() + 
  xlab(NULL) + 
  ylab("Count") + 
  labs(title = "PAD2")  + 
  ylim(0,60) 

  
  ggplot(df4, aes(fill = cit, x = type)) + 
  geom_bar(stat = "count") + 
  theme_prism() + 
  xlab(NULL) + 
  ylab("Count") + 
  labs(title = "PAD4") + 
  ylim(0,60) 

#Fig B. Violin plot of distance to cit for created and destroyed regions.----
  
  ggplot(df1, aes(fill = cit, x = type)) + 
    geom_bar(stat = "count") + 
    theme_prism() + 
    xlab(NULL) + 
    ylab("Count") + 
    labs(title = "Cits in DiffEx Regions") + 
    #ylim(0,60) + 
    facet_grid(~ cond)
  
  
  ggplot(brdf, aes(x = type, y = diff)) + 
    geom_violin() + 
    theme_prism() + 
    ylab("Actual - Exp dist") + 
    facet_grid(~cond)
  
#  save(brdf, file = "./objects/regions8Apr.rda")
#Fig C. BoxPlot of local variances. 

#Save data
 region_data <- brdf
 save(region_data, "./objects/region_data.rda")