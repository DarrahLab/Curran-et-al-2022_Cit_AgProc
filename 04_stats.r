#Some statistics. 
#A) For dist to cit, show convincingly the distributions are not normal. 
# > Shapiro Wilk Test
# > QQ Plot. 
#B)

library(stats)
library(gridExtra)
library(ggplot2)
library(diptest)
setwd("C:/Users/aag73/Desktop/ZZZ_protmap")
load("./objects/13MarStuff.rda")
load("./objects/regions8Apr.rda")
load("./objects/varsdf_8Apr.rda")
load("./objects/protmap_dat.rda")
load("./objects/cit-sites-all-10.rda")
load("./objects/diffExdat6Apr.rda")
region_data <- brdf
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
variances_df <- big_df

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
t <- c(1,2,3,5)
diffEx <- diffEx[t]
protmap_dat <- protmap_dat[t]

#load sequences. -----
seq.ra33<- read.csv("./fibrinogen_protmap/seq.ra33.csv", header = FALSE)$V1
seq.alpha <- read.csv("./fibrinogen_protmap/seq.fiba.csv", header = FALSE)$V1
seq.beta <- read.csv("./fibrinogen_protmap/seq.fibb.csv", header = FALSE)$V1
seq.vim<- read.csv("./fibrinogen_protmap/seq.vim.csv", header = FALSE)$V1
seq.pad4<- read.csv("./fibrinogen_protmap/seq.pad4.csv", header = FALSE)$V1
seq.gamma <-  read.csv("./fibrinogen_protmap/seq.fibg.csv", header = FALSE)$V1

seq_list <- list(seq.ra33, seq.alpha,seq.beta , seq.vim,seq.pad4 )
names(seq_list) <- c("ra33", "fib.a", "fib.b", "vim", "pad4")
rm(seq.alpha,seq.beta,seq.pad4,seq.ra33,seq.vim)

for(i in 1:length(seq_list)){
  seq_list[[i]] <- strsplit(seq_list[[i]], "")[[1]]
}

#Ascertain protein coverage - 'deep' reads####
#Find sites in pad2 and pad4 conditions abd > cutoff (1e6) We say these 
#are read at enough depth to call a citrulline if it twere citrullinated. 
Ags <- names(protmap_dat)
deep_reads <- vector(mode = "list", length = 4); names(deep_reads) <- Ags
cutoff <- 1e6
seq_depth <- matrix(data = NA, nrow = 4, ncol = 2)
colnames(seq_depth) <- c("pad2", "pad4")
rownames(seq_depth) <- Ags
for(a in 1:length(Ags)){
  tmp <- list(pad2 = NA, pad4 = NA)
  
  p2.deep <- which(protmap_dat[[a]][[2]] >= cutoff)
  tmp[[1]] <- p2.deep
  seq_depth[a,1] <- length(p2.deep) /  length(seq_list[[a]])
  
  if(a < 6){
    p4.deep <- which(protmap_dat[[a]][[3]] >= cutoff)
    tmp[[2]] <- p4.deep
    seq_depth[a,2] <- length(p4.deep) /  length(seq_list[[a]])
  }
  
  deep_reads[[a]] <- tmp
}
rm(tmp, p2.deep, p4.deep)
#D to Cit SW Test and QQ Plot----

#Split up region data into pad2/4, created/destroyed as indicated in plots. 
c2idx <- which(region_data$cond == "pad2" & region_data$type == "Created")
d2idx <- which(region_data$cond == "pad2" & region_data$type == "Destroyed")
c4idx <- which(region_data$cond == "pad4" & region_data$type == "Created")
d4idx <- which(region_data$cond == "pad4" & region_data$type == "Destroyed")

nms <- c("Pad2 Created", "Pad2 Destroyed", "Pad4 Created", "Pad4 Destroyed")
idx <- list(c2idx, d2idx, c4idx, d4idx); names(idx) <- nms
swStat <- vector(mode = "list", length = 4); names(swStat) <- nms
dipStat <- vector(mode = "list", length = 4); names(swStat) <- nms
qq <- vector(mode = "list", length = 4); names(qq) <- nms
for(i in 1:4){
  d <- region_data$diff[idx[[i]]]
  swStat[[i]] <- shapiro.test(d)$p.value
  dipStat[[i]] <- dip(d)
  qq[[i]] <- qqnorm(d)
  #qqnorm(d[which(d < quantile(d,0.85))])
}

sp <- vector(mode = "list", length = 4)

for(p in 1:4){
  d <- qq[[p]]
  df <- data.frame(x = d$x, y = d$y)
  slope <- (max(df$y) - min(df$y)) / (max(df$x) - min(df$x))
  int <- max(df$y) - slope * max(df$x)
  
  
  g <- ggplot(data = df, mapping = aes(x = x, y = y)) + 
       geom_point() + 
       geom_abline(slope = slope, intercept = int , color = "dark grey", linetype = "dashed") + 
      theme_prism() + 
      xlab("Theoretical Quantiles") + 
    ylab("Sample Quantiles") + 
    labs(subtitle = nms[[p]] )
  
  sp[[p]] <- g
}

grid.arrange(sp[[1]], sp[[2]], sp[[3]], sp[[4]], ncol = 2)



d <- qqnorm(region_data$diff)
df <- data.frame(x = d$x, y = d$y)
slope <- (max(df$y) - min(df$y)) / (max(df$x) - min(df$x))
int <- max(df$y) - slope * max(df$x)


g <- ggplot(data = df, mapping = aes(x = x, y = y)) + 
  geom_point() + 
  geom_abline(slope = slope, intercept = int , color = "dark grey", linetype = "dashed") + 
  theme_prism() + 
  xlab("Theoretical Quantiles") + 
  ylab("Sample Quantiles") + 
  labs(subtitle = nms[[p]] )

#Variance sig testing. -------
variances_df <- variances_df[-which(variances_df$Ag == "pad4"),]
#Run bootstrap, get pval. 
N <- 500
Ags <- unique(variances_df$Ag)
nAg <- length(Ags)
pvals <- vector(mode = "double", length = nAg); names(pvals) <- Ags
ratios <- vector(mode = "double", length = nAg); names(ratios) <- Ags

#Same but with medians
results  <- data.frame(group = NA, median = NA, med_exp = NA,med_norm = NA, stdev = NA,
                       pval = NA, sortby = NA )
#Broken up by Antigen 
for (a in 1:length(Ags)){
  dat <- variances_df$Val[which(variances_df$Ag == Ags[a])]
  cit <- variances_df$Val[which(variances_df$Ag == Ags[a] & variances_df$Class == "Cit")]
  bootidx <- rcomb(length(dat), length(cit), N)
  
  meds <- vector(mode = "double", length = N)
  for(b in 1:N){
    meds[b] <- median(dat[bootidx[[b]]])
  }
  
  tmp <- data.frame(group = paste0(Ags[a], "_all"),
                    median = median(cit), 
                    med_exp = mean(meds), 
                    med_norm = median(cit) / mean(meds),
                    pval = sum(meds >= median(cit)) / N, 
                    stdev = sd(meds) / mean(meds),
                    sortby = "ag")
  results<- rbind(results,tmp)
  
}

#Broken by antigen and pad2/pad4
conds <- c("pad2", "pad4")
for (a in 1:length(Ags)){
  for(c in 1:2){
    dat <- variances_df$Val[which(variances_df$Ag == Ags[a])]
    cit <- variances_df$Val[which(variances_df$Ag == Ags[a] & variances_df$Class == "Cit" 
                                  &variances_df$Cond == conds[c])]
    bootidx <- rcomb(length(dat), length(cit), N)
    
    meds <- vector(mode = "double", length = N)
    for(b in 1:N){
      meds[b] <- median(dat[bootidx[[b]]])
    }
    
    tmp <- data.frame(group = paste0(Ags[a], "_",conds[c]),
                      median = median(cit), 
                      med_exp = mean(meds), 
                      med_norm = median(cit) / mean(meds),
                      pval = sum(meds >= median(cit)) / N,
                      stdev = sd(meds) / mean(meds),
                      sortby = "ag_pad")
    results <- rbind(results,tmp)
    
  }
}

#Broken up by condition
conds <- c("pad2", "pad4")

for(c in 1:2){
  dat <- variances_df$Val[which(variances_df$Cond == conds[c])]
  cit <- variances_df$Val[which( variances_df$Class == "Cit" &variances_df$Cond == conds[c])]
  bootidx <- rcomb(length(dat), length(cit), N)
  
  meds <- vector(mode = "double", length = N)
  for(b in 1:N){
    meds[b] <- median(dat[bootidx[[b]]])
  }
  
  tmp <- data.frame(group = paste0(conds[c], "_all"),
                    median = median(cit), 
                    med_exp = mean(meds), 
                    med_norm = median(cit) / mean(meds),
                    pval = sum(meds >= median(cit)) / N, 
                    stdev = sd(meds) / mean(meds),
                    sortby = "pad")
  results <- rbind(results,tmp)
  
}
results <- results[-1,]

median_results <- results

#Plotting var sig tests------
tmp <- median_results[which(results$sortby == "ag_pad"),]
tmp$pad <- NA
tmp$pad[grep("pad2", tmp$group)] <- "pad2"
tmp$pad[grep("pad4", tmp$group)] <- "pad4"
ggplot(tmp, mapping = aes(x = group, y = med_norm, fill = pad)) + 
  geom_col() + 
  theme_prism() + 
  geom_hline(yintercept = 1, color = "red") + 
  theme(axis.text.x = element_text(angle = 45))



hist(variances_df$Val[which(variances_df$Class == "Not")])

d <- variances_df$Val[which(variances_df$Class == "Not")]
e <-  variances_df$Val[which(variances_df$Class == "Cit")]
ggplot(data.frame(x = d), aes(x = x, y = ..density..)) +
  geom_density() +
  geom_density(data = data.frame(x = e), color = "red", alpha = 0.5)

qqplot(d,e)


d2 <- d[which(d <= quantile(d, 0.4))]
e2 <- e[which(e <= quantile(e, 0.4))]
qqplot(d2,e2)


#Chi square - distribution of cits in Cr, De, NC regions ----
region_data <- split(region_data,region_data$ag)
region_data <- region_data[c(4,1,2,5)]

conds <- c("pad2", "pad4")


pad2 <- vector(mode = "integer", length = 6)
pad4 <- vector(mode = "integer", length = 6)
names(pad2) <- c("CrRes", "DeRes", "NCRes", "CrCit", "DeCit", "NCCit")
names(pad4) <- names(pad2)
for(a in 1:4){
  for(p in 1:2){
    totRes <- length(deep_reads[[a]][[p]])
    nCit <- length(cit_sites_all[[p]][[a]])
    d <- region_data[[a]][which(region_data[[a]]$cond == conds[p]),]
    Cr <- d[which(d$type == "Created"),]
    De <- d[which(d$type == "Destroyed"),]
    CrRes <- sum(Cr$end - Cr$start + 1)
    DeRes <- sum(De$end - De$start + 1)
    NCRes <- totRes - CrRes - DeRes
    CrCit <- sum(Cr$dist == 0)
    DeCit <- sum(De$dist == 0)
    NCCit <- nCit - CrCit - DeCit
    
    if(p == 1){
      pad2 <- pad2 + c(CrRes,DeRes,NCRes,CrCit,DeCit,NCCit)
    }else{
      pad4 <- pad4 + c(CrRes,DeRes,NCRes,CrCit,DeCit,NCCit)
    }
  }
}


pad2dat <- matrix(nrow = 2, ncol = 3)
colnames(pad2dat) <- c("Cr", "De", "NC")
rownames(pad2dat) <- c("Not", "Cit")

pad2dat[1,] <- c(pad2[1] - pad2[4], 
                 pad2[2] - pad2[5], 
                 pad2[3] - pad2[6])
pad2dat[2,] <- c( pad2[4], 
                  pad2[5], 
                  pad2[6])




pad4dat <- matrix(nrow = 2, ncol = 3)
colnames(pad4dat) <- c("Cr", "De", "NC")
rownames(pad4dat) <- c("Not", "Cit")

pad4dat[1,] <- c(pad4[1] - pad4[4], 
                 pad4[2] - pad4[5], 
                 pad4[3] - pad4[6])
pad4dat[2,] <- c( pad4[4], 
                  pad4[5], 
                  pad4[6])



#RNG for bootstrap
bootagain <- function(dat){
  N = 1000
  nRes <- colSums(dat)
  nCit <- rowSums(dat)[2]
  allRes <- sum(nRes)
  boots <- rcomb(allRes, nCit, N)
  d_exp <- matrix(nrow = N, ncol = 3)
  for(b in 1:N){
    cits <- as.double(boots[[b]])
    nCr <- sum(cits <= nRes[1])
    nDe <- sum(cits > nRes[1] & cits <= (nRes[1] + nRes[2]))
    nNC <- sum(cits > (nRes[1] + nRes[2]))
    
    d_exp[b,] <- c(nCr, nDe, nNC)
  }
  
  pCr <- sum(d_exp[,1] >= dat[2,1]) / N
  pDe <- sum(d_exp[,2] >= dat[2,2]) / N
  pNC <- sum(d_exp[,3] >= dat[2,3]) / N
  
  out <- c(pCr, pDe, pNC)
  names(out) <- c("Cr", "De", "NC")
  return(out)
}


bootagain(pad2dat)
bootagain(pad4dat)



# 
# nCit2 <- sum(pad2[4:6])
# CrRes <- pad2[1]; DeRes <- pad2[2]; NCRes <- pad2[3]
# CrE2 <- CrRes / (CrRes + DeRes + NCRes) * nCit2
# DeE2 <- DeRes/ (CrRes + DeRes + NCRes) * nCit2
# NCE2 <- NCRes / (CrRes + DeRes + NCRes) * nCit2
# 
# pad2dat <- matrix(nrow = 3, ncol = 2)
# pad2dat[,2] <- pad2[4:6]
# pad2dat[,1] <- c(CrE2, DeE2, NCE2)
# colnames(pad2dat) <- c("Exp", "Actual")
# rownames(pad2dat) <- c("Cr", "De", "NC")
# 
# 
# nCit4 <- sum(pad4[4:6])
# CrRes <- pad4[1]; DeRes <- pad4[2]; NCRes <- pad4[3]
# CrE4 <- CrRes / (CrRes + DeRes + NCRes) * nCit4
# DeE4 <- DeRes/ (CrRes + DeRes + NCRes) * nCit4
# NCE4 <- NCRes / (CrRes + DeRes + NCRes) * nCit4
# 
# pad4dat <- matrix(nrow = 3, ncol = 2)
# pad4dat[,2] <- pad4[4:6]
# pad4dat[,1] <- c(CrE4, DeE4, NCE4)
# colnames(pad4dat) <- c("Exp", "Actual")
# rownames(pad4dat) <- c("Cr", "De", "NC")
# 
# library(stats)
# chisq.test(x = pad2dat[,2], y = pad2dat[,1])
# chisq.test(x = pad4dat[,2], y = pad4dat[,1])