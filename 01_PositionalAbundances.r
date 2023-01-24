##From mass spec data, obtain positional abundance by residue for each antigen. 
## 20 February 2022 AAG

rm(list = ls())
###LIBRARIES/WORKDIR####
setwd("/media/aag7319/WDRed/Curran-et-al-2022_Cit_AgProc")
library(readxl)
library(strex)
library(seqinr)
source('./functions/fn.load.r')
fn.load()
###LOAD,FORMAT DATA####
run1 <- read.all.sheets("./SourceData/CatBSH_digested_RA33_results_042921.xlsx")
run2 <- read.all.sheets(paste0("./SourceData/CatBSH_digested_fibrinogen",
                               "_vimentin_PAD4_results_072821.xlsx"))
run3 <- read.all.sheets("./SourceData/05_2021_12_28.xlsx")

v.interest <- c("Sequence", "Modifications", "Master Protein Accessions", 
                "Positions in Master Proteins")

#Extract from run1
s.idx <- c(15:17, 24:26, 18:20, 27:29, 21:23, 30:32)
d.idx <- c(match(v.interest, colnames(run1[[2]])), s.idx )

dat1 <- run1[[2]][,d.idx ]
dat1 <- dat1[(dat1$`Master Protein Accessions` == "P22626"), ] #subset to RA33
dat1$Protein <- "RA33"; dat1$file <- "run1"

#Extract from run2
dat2 <- list(fibrinogen = NA, vimentin = NA, pad4 = NA)
#fibrinogen
s.idx <- c(19,22,25,28,20,23,26,29,21,24,27,30)
d.idx <- c(match(v.interest, colnames(run2[[2]])), s.idx )
dat2[[1]] <- run2[[2]][,d.idx]
p.idx <- c( grep("P02671", dat2[[1]]$`Master Protein Accessions` ), 
            grep("P02675", dat2[[1]]$`Master Protein Accessions` ),
            grep("P02679", dat2[[1]]$`Master Protein Accessions` ))

dat2[[1]] <- dat2[[1]][p.idx, ]
#vimentin 
s.idx <- c(19,22,25,28,20,23,26,29,21,24,27,30) - 1
d.idx <- c(match(v.interest, colnames(run2[[4]])), s.idx )
dat2[[2]] <- run2[[4]][,d.idx]
p.idx <- grep("P08670", dat2[[2]]$`Master Protein Accessions` )
dat2[[2]] <- dat2[[2]][p.idx, ]
#PAD4
s.idx <-c(19,21,23,25, 20,22,24,26)
d.idx <- c(match(v.interest, colnames(run2[[6]])), s.idx )
dat2[[3]] <- run2[[6]][,d.idx]
p.idx <- grep("Q9UM07", dat2[[3]]$`Master Protein Accessions` )
dat2[[3]] <- dat2[[3]][p.idx, ]

#Extract from run3
dat3 <- list(fibrinogen = NA, vimentin = NA, pad4 = NA)
colnames(run3[[2]]) <- run3[[2]][1,]; # run3[[2]][-1,]
colnames(run3[[2]])[19:26] <- colnames(run3[[1]])[18:25]
#fibrinogen
s.idx <- c(19:21)
d.idx <- c(match(v.interest, colnames(run3[[2]])), s.idx )
p.idx <- c( grep("P02671", run3[[2]]$`Master Protein Accessions` ), 
            grep("P02675", run3[[2]]$`Master Protein Accessions` ),
            grep("P02679", run3[[2]]$`Master Protein Accessions` ))
dat3[[1]]<- run3[[2]][p.idx,d.idx]
#vimentin
s.idx <- c(22:24)
d.idx <- c(match(v.interest, colnames(run3[[2]])), s.idx )
p.idx <- grep("P08670", run3[[2]]$`Master Protein Accessions` )
dat3[[2]] <- run3[[2]][p.idx,d.idx]
#pad4
s.idx <- c(25:26)
d.idx <- c(match(v.interest, colnames(run3[[2]])), s.idx )
p.idx <- grep("Q9UM07", run3[[2]]$`Master Protein Accessions` )
dat3[[3]] <- run3[[2]][p.idx,d.idx]
###CLEAN DATA####
#obtain master sequence starts and ends. 
tmp <- bookEnds(dat1$`Positions in Master Proteins`)
dat1$Start <- tmp[,1]; dat1$End <- tmp[,2]; rm(tmp)
acc <-  c(NULL, "P08670", "Q9UM07")
for (i in 1:3){
  
  if(i == 3){
    tmp2 <- bookEnds(dat2[[i]]$`Positions in Master Proteins`, "Q9UM07")
    tmp3 <- bookEnds(dat3[[i]]$`Positions in Master Proteins`,  "Q9UM07")
  }else if(i == 2){
    tmp2 <- bookEnds(dat2[[i]]$`Positions in Master Proteins`, "P08670")
    tmp3 <- bookEnds(dat3[[i]]$`Positions in Master Proteins`,  "P08670")
  }else {
    tmp2 <- bookEnds(dat2[[i]]$`Positions in Master Proteins`)
    tmp3 <- bookEnds(dat3[[i]]$`Positions in Master Proteins`) 
  }
  
  dat2[[i]]$Start <- tmp2[,1];  dat2[[i]]$End <- tmp2[,2]
  dat3[[i]]$Start <- tmp3[,1];  dat3[[i]]$End <- tmp3[,2]
  rm(tmp2,tmp3)
}

#label file of origin prior to conacatenation.
dat1 <- only.cit(dat1)
for (i in 1:3){
  dat2[[i]]$file = "run2"
  dat3[[i]]$file = "run3"
}

###SUBSET BY CONDITION, GET CIT PCT####
#RA33, the easy one. 
ra33.pcts <- list(native = positional.abundance(dat1, c(5:10), "ra33"), 
                  pad2 = positional.abundance(dat1, c(11:16), "ra33") ,
                  pad4 = positional.abundance(dat1, c(17:22), "ra33") )
#Fibrinogen
fib.acc <- c("P02671","P02675", "P02679") #alpha, beta, gamma accession #s
tmp <- list(native = NA, pad2 = NA, pad4 = NA)
fib.pcts <- list(alpha = tmp, beta = tmp, gamma = tmp)
fib.run2.pcts <- list(alpha = tmp, beta = tmp, gamma = tmp)
fib.run3.pcts <- list(alpha = tmp, beta = tmp, gamma = tmp)
#dat2: 4 replicates. dat3: no replicates. 
#dat.idx is interpersed locations of nat.dat2, nat.dat3, pad2.dat2, pad2.dat3 etc. 
dat.idx <- list(c(5:8), c(5), c(9:12), c(6), c(13:16), c(7))
for(c in 1:3){ #for each condition 
  for (s in 1:3){ #for each subunit 
    #subset data from dat2
    d2.idx <- which(dat2$fibrinogen$`Master Protein Accessions` == fib.acc[s])
    d2 <- dat2[[1]][d2.idx,]
    d2.pct <- positional.abundance(d2, dat.idx[[2*c - 1]], fib.acc[s])
  
    #subset data from dat3
    d3.idx <- which(dat3$fibrinogen$`Master Protein Accessions` == fib.acc[s])
    d3 <- dat3[[1]][d3.idx,]
    d3.pct <- positional.abundance(d2, dat.idx[[2*c]], fib.acc[s] )
    
    fib.run2.pcts[[s]][[c]]<- d2.pct 
    fib.run3.pcts[[s]][[c]]<- d3.pct
  }}

#Vimentin
vim.pcts <- list(native = NA, pad2 = NA, pad4 = NA)
vim.run2.pcts <- list(native = NA, pad2 = NA, pad4 = NA)
vim.run3.pcts <- list(native = NA, pad2 = NA, pad4 = NA)
dat.idx <- list(c(5:8), c(5), c(9:12), c(6), c(13:16), c(7))

for(c in 1:3){
  d2.pct <- positional.abundance(dat2[[2]], dat.idx[[2*c - 1]], "vim")
  d3.pct <- positional.abundance(dat3[[2]], dat.idx[[2 * c]], "vim")
  vim.run2.pcts[[c]] <- d2.pct
  vim.run3.pcts[[c]] <- d3.pct
  vim.pcts[[c]] <- (d2.pct + d3.pct) / 2
}

#pad4 
pad4.pcts <- list(native = NA, cit = NA)
pad4.run2.pcts <- list(native = NA, cit = NA)
pad4.run3.pcts <- list(native = NA, cit = NA)
dat.idx <- list(c(5:8), c(5), c(9:12), c(6) )

for(c in 1:2){
  d2.pct <- positional.abundance(dat2[[3]], dat.idx[[2*c-1]], "pad4")
  d3.pct <- positional.abundance(dat3[[3]], dat.idx[[2*c]], "pad4")
  pad4.pcts[[c]] <- (d2.pct + d3.pct) / 2
  pad4.run2.pcts[[c]] <- d2.pct 
  pad4.run3.pcts[[c]] <- d3.pct
}


data <- list(fiba = fib.pcts[[1]], fibb = fib.pcts[[2]], 
             fibg = fib.pcts[[3]], vim = vim.pcts, ra33 = ra33.pcts,
             pad4 = pad4.pcts)

data.run2 <- list(fiba = fib.run2.pcts[[1]], fibb = fib.run2.pcts[[2]], 
                  fibg = fib.run2.pcts[[3]], vim = vim.run2.pcts,
                  pad4 = pad4.run2.pcts)

data.run3 <- list(fiba = fib.run3.pcts[[1]], fibb = fib.run3.pcts[[2]], 
                  fibg = fib.run3.pcts[[3]], vim = vim.run3.pcts,
                  pad4 = pad4.run3.pcts)

#Group data by native or pad condition, NOT averaging runs but just throwing them 
#all together. used in scatter by condition section. 
data.native <- as.data.frame(data[[5]][[1]])
data.native$Ag <- "ra33"
data.native$residue <- c(1:353)
data.native$run <- "Run1"

data.pad2 <- as.data.frame(data[[5]][[2]])
data.pad2$Ag <- "ra33"
data.pad2$residue <- c(1:353)
data.pad2$run <- "Run1"
data.pad2$cond <- "pad2"

data.pad4 <- as.data.frame(data[[5]][[3]])
data.pad4$Ag <- "ra33"
data.pad4$residue <- c(1:353)
data.pad4$run <- "Run1"
data.pad4$cond <- "pad4"
aglen <- c(866,491,453,466,663)
ags <- c("fib.a", "fib.b", "fib.g", "vim", "pad4")
for (i in 1:5){
  x <- as.data.frame.matrix(data.run2[[i]][[1]] )
  y <- as.data.frame.matrix(data.run3[[i]][[1]] )
  x$residue <-c(1:aglen[i]) ; x$Ag <- ags[i]; x$run <- "Run2"
  y$residue <- c(1:aglen[i]); y$Ag <- ags[i]; y$run <- "Run3"
  data.native <- rbind(data.native, x, y)
  
  x <- as.data.frame.matrix(data.run2[[i]][[2]] )
  y <- as.data.frame.matrix(data.run3[[i]][[2]] )
  x$residue <-c(1:aglen[i]) ; x$Ag <- ags[i]; x$run <- "Run2"
  y$residue <- c(1:aglen[i]); y$Ag <- ags[i]; y$run <- "Run3"
  x$cond <- "pad2"
  y$cond <- "pad2"
  data.pad2 <- rbind(data.pad2, x, y)
  
  if (i < 5){
    x <- as.data.frame.matrix(data.run2[[i]][[3]] )
    y <- as.data.frame.matrix(data.run3[[i]][[3]] )
    x$residue <-c(1:aglen[i]) ; x$Ag <- ags[i]; x$run <- "Run2"
    y$residue <- c(1:aglen[i]); y$Ag <- ags[i]; y$run <- "Run3"
    x$cond <- "pad4"
    y$cond <- "pad4"
    data.pad4<- rbind(data.pad4, x, y)
  }
}

readme_protmap_dat <- "Protmap_dat contains UNadjusted abundances and %cit mass spec data from 07/28/21 (all Ag except ra33) and 04/29/21 (ra33). "
data.run2[[6]] <- ra33.pcts

protmap_dat <- data.run2[c(6,1:5)];  names(protmap_dat)[1] <- "ra33"
rm(d2,d2.pct,d3,d3.pct,dat.idx,dat1,dat2,dat3,run1,run2,run3,tmp, 
   acc,c,d.idx,d2.idx,d3.idx,fib.acc,i,p.idx,v.interest)
save(protmap_dat, file = "./objects/protmap_dat.rda")

