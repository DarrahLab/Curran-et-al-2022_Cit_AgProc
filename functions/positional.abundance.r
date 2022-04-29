positional.abundance <- function(df.in, abd.idx, id){
  
  st.idx <- grep("Start", colnames(df.in))[1]
  sp.idx <- grep("End", colnames(df.in))[1]
  if (is.null(st.idx) || is.null(sp.idx)){
    print("Error: Improperly formatted input column names.")
    return(NULL)
  }
  
  # Accession Numbers for Reference: 
    hmgcr <- c("P04035", "hmgcr")
    fib.a <- c("P02671", "fib.a")
    fib.b <- c("P02675", "fib.b")
    fib.g <- c("P02679", "fib.g")
    ra33 <- c("P22626", "ra33")
    vim <- c("P08670", "vim")
    pad4 <- c("Q9UM07", "pad4")
  
  #Determine total length of protein
  if (id %in% hmgcr){ len <- 888} else
  if (id %in% fib.a){ len <- 866} else 
  if (id %in% fib.b){ len <- 491} else 
  if (id %in% fib.g){ len <- 453} else
  if (id %in% ra33){ len <- 353} else 
  if (id %in% vim ){len <- 466} else 
  if (id %in% pad4){len <- 663} else{
    print("Error: Protein input not recognized.")
    return(NULL)
  }
    
  master <- vector(mode = "double", length = len)
  
  for (i in 1:dim(df.in)[1]){
    tot.abd <- sum(as.numeric(df.in[i,abd.idx]), na.rm = TRUE)
    #Add average abundance. 
    master[df.in[i,st.idx]:df.in[i,sp.idx]] = tot.abd/length(abd.idx) + 
                                                master[df.in[i,st.idx]:df.in[i,sp.idx]] 
  }
  
  master[is.na(master)] = 0
  print(paste("Abundance count successfully returning...",id,sep = ""))
  return(master)
    
}
positional.abd.cit <- function(df.in, abd.idx, id){
  
  st.idx <- grep("Start", colnames(df.in))[1]
  sp.idx <- grep("End", colnames(df.in))[1]
  if (is.null(st.idx) || is.null(sp.idx)){
    print("Error: Improperly formatted input column names.")
    return(NULL)
  }
  
  # Accession Numbers for Reference: 
  hmgcr <- c("P04035", "hmgcr")
  fib.a <- c("P02671", "fib.a")
  fib.b <- c("P02675", "fib.b")
  fib.g <- c("P02679", "fib.g")
  ra33 <- c("P22626", "ra33")
  vim <- c("P08670", "vim")
  pad4 <- c("Q9UM07", "pad4")
  
  #Determine total length of protein
  seq <- NA
  if (id %in% hmgcr){ len <- 888} else
    if (id %in% fib.a){ 
     len <- 866; 
     seq <- read.csv("./SourceData/seq.fiba.csv", header = FALSE)$V1 
    } else if (id %in% fib.b){ 
      len <- 491
      seq<- read.csv("./SourceData/seq.fibb.csv", header = FALSE)$V1
    } else if (id %in% fib.g){ 
      len <- 453
      seq <- read.csv ("./SourceData/seq.fibg.csv", header = FALSE)$V1
    } else if (id %in% ra33){ 
      len <- 353
      seq <- read.csv("./SourceData/seq.ra33.csv", header = FALSE)$V1
    } else if (id %in% vim ){
      len <- 466
      seq <- read.csv ("./SourceData/seq.vim.csv", header = FALSE)$V1
    } else if (id %in% pad4){
      len <- 663
      seq <- read.csv("./SourceData/seq.pad4.csv", header = FALSE)$V1
      } else{
      print("Error: Protein input not recognized.")
      return(NULL)
      }
  
  seq_char <-  strsplit(seq, split = character(0))[[1]]
  master <- vector(mode = "double", length = len)
  
  for (i in 1:dim(df.in)[1]){
    is.cit <- length(grep("Deamidated", df.in$Modifications[i])) > 0
    if(!is.cit){next}else{
    
    cit.sites <- find_deamidation(df.in[i,])
    
    if(!all(seq_char[cit.sites] == "R")){
      print("Warning: Non-Arginine Deamidations Detected! Filtering...")
      print(paste0("i = ", i))
      cit.sites <- cit.sites[which(cit.sites == "R")]
    }
    
    tot.abd <- sum(as.numeric(df.in[i,abd.idx]), na.rm = TRUE)
    #Add average abundance. 
    master[cit.sites] = tot.abd/length(abd.idx) + 
      master[cit.sites] 
  }}
  
  master[is.na(master)] = 0
  print(paste("Cit count successfully returning...",id,sep = ""))
  return(master)
  
}
positional.cit.pct <- function(df.in, abd.idx, id){
  pos.abd <- positional.abundance(df.in, abd.idx, id)
  pos.abd.cit <- positional.abd.cit(df.in,abd.idx,id)
  cit.pct <- pos.abd.cit
  cit.idx <- which(pos.abd.cit != 0)
  cit.pct[cit.idx] <- pos.abd.cit[cit.idx] / pos.abd[cit.idx] * 100
  
  #pos.abd <- pos.abd / sum(pos.abd) * 100 #present adjusted abundance
  out <- cbind(cit.pct, pos.abd)
  colnames(out) <- c("Percent Citrullinated", "Total Abd")
  return(out)
}
plot.cit.pct <- function(dat, cond, specifier){
  #Input is an r x 2 matrix, with column 1 names "Percent Citrullinated"
  #and column 2 "Total Abd"
  #Plot them both to visualize citrullination locations and magnitude
  #in your sample 
  
  if(missing(cond)){
    cond <- NULL
  }
  
  if(missing(specifier)){
    specifier<- NULL
  }else{
    specifier <- paste0(": ",specifier )
  }
  #Determine protein (for plot title)
  l <- dim(dat)[1]
  if(l == 866){p <- "fib.a"}else
    if(l == 491){p <- "fib.b"}else
      if(l == 453){p <- "fib.g"}else
        if(l == 353){p <- "ra33"}else
          if(l == 466){p <- "vim"}else
            if(l == 663){p <- "pad4"}else
            {p <- "unknown protein"}
  df <- data.frame(dat)
  df$pos <- c(1:dim(df)[1])
  
  
  
  scale.factor <- (max(df$Percent.Citrullinated) / 2) / 
    max(df$Total.Abd)
  df$Total.Abd <- df$Total.Abd * scale.factor
  
  df1 <- df[,-2] #cit pcts
  df1 <- df1[which(df1[,1] != 0),]
  df2 <- df[,-1] #tot abd
  colnames(df1) <- c("y", "pos")
  colnames(df2) <- c("y", "pos")
  
  p <- ggplot(df1, aes(x = pos, y = y)) + geom_point(size = 1, fill = "black", shape = 15) + 
    geom_col(width = 1, fill = "black") + theme_minimal() + xlab("AA Position") +
    geom_area(data = df2, alpha = 0.3, color = mako(1, begin = 0.5), 
              fill = mako(1,begin = 0.5), size = 0.2) + ylab("% Citrullinated") + 
    labs(title = paste0(p, " ", cond,specifier )) + ylim(0,max(df$Percent.Citrullinated))
  return(p)
}
plot.cit.hist <- function(dat, cond, specifier){
  #Input is an r x 2 matrix, with column 1 names "Percent Citrullinated"
  #and column 2 "Total Abd"
  #specifier is a string to add to title. 
  #Plot them both to visualize citrullination locations and magnitude
  #in your sample 
  
  if(missing(cond)){
    cond <- NULL
  }
  
  if(missing(specifier)){
     specifier<- NULL
  }else{
    specifier <- paste0(": ",specifier )
  }
  #Determine protein (for plot title)
  l <- dim(dat)[1]
  if(l == 866){p <- "fib.a"}else
    if(l == 491){p <- "fib.b"}else
      if(l == 453){p <- "fib.g"}else
        if(l == 353){p <- "ra33"}else
          if(l == 466){p <- "vim"}else
            if(l == 663){p <- "pad4"}else
            {p <- "unknown protein"}
 
   df <- data.frame(dat)
  df$pos <- c(1:dim(df)[1])
  
  df1 <- df[,-2] #cit pcts
  df1 <- df1[which(df1[,1] != 0),]
  df2 <- df[,-1] #tot abd
  colnames(df1) <- c("y", "pos")
  colnames(df2) <- c("y", "pos")
  
  p <- ggplot(df1, aes(x = y)) + ylab("# Residues") + xlab("% Citrullinated") + 
    labs(title = paste0(p, " ", cond, specifier )) + geom_histogram() + theme_minimal()
  return(p)
}
plot.cit.scatter <- function(data, specifier){
  #Input is an list of 2 or 3 cond: native, pad2, pad4 for a given ag.  
  # for each item in list, 3 r x 2 matrices with column 1 names "Percent Citrullinated"
  #and column 2 "Total Abd"
  #specifier is a string to add to title. 
  #Plot them both to visualize citrullination locations and magnitude
  #in your sample 
  
  if (length(data) == 3){
    LEN = "THREE"
    cond = c("native", "pad2", "pad4")
  }else if (length(data) == 2){
    LEN = "TWO"
    cond = c("native", "cit")
  }else{
    print("Error: Unexpected input list length!")
    return(NULL)
  }
  
  if(missing(specifier)){
    specifier<- NULL
  }else{
    specifier <- paste0(": ",specifier )
  }
  #Determine protein (for plot title)
  l <- dim(data[[1]])[1]
  if(l == 866){p <- "fib.a"}else
    if(l == 491){p <- "fib.b"}else
      if(l == 453){p <- "fib.g"}else
        if(l == 353){p <- "ra33"}else
          if(l == 466){p <- "vim"}else
            if(l == 663){p <- "pad4"}else
            {p <- "unknown protein"}
  
  dfs <- vector(mode = "list", length = length(data))
  for(i in 1:length(data)){
    dat <- data[[i]]
    zero.idx <- which(dat[,1] == 0)
    dat <- dat[-zero.idx,]
    
    df <- try(as.data.frame.matrix(dat), silent = TRUE)
    if(typeof(df) == "character"){ df <- data.frame(dat)}
    dfs[[i]] <- df
  }
  
  if (LEN == "THREE"){
    p0 <- ggplot(dfs[[1]], aes(x = `Total Abd`, y = `Percent Citrullinated`)) + 
      ylab("% Citrullinated") + xlab("log(Abd)") + 
      labs(title = paste0(p, " ", specifier )) + geom_point( color = "black", shape = 1, fill = "white") + 
      geom_point(data = dfs[[2]], color = "red") + geom_point(data = dfs[[3]], color = "black") + 
      theme_classic() + theme(panel.grid.major = element_line(colour = "light grey")) +  scale_x_log10() 
  }else if (LEN == "TWO"){
    p0 <- ggplot(dfs[[1]], aes(x = `Total Abd`, y = `Percent Citrullinated`)) + 
      ylab("% Citrullinated") + xlab("log(Abd)") + 
      labs(title = paste0(p, " ", specifier )) + geom_point( color = "black", shape = 1, fill = "white") + 
      geom_point(data = dfs[[2]], color = "black")  + 
      theme_classic() + theme(panel.grid.major = element_line(colour = "light grey")) +  scale_x_log10() 
  }

  
  return(p0)
}
plot.cit.scatter2 <- function(data.run2,data.run3){
  #Input is an list of 2 or 3 cond: native, pad2, pad4 for a given ag.  
  # for each item in list, 3 r x 2 matrices with column 1 names "Percent Citrullinated"
  #and column 2 "Total Abd"
  #specifier is a string to add to title. 
  #Plot them both to visualize citrullination locations and magnitude
  #in your sample 
  

  #Determine protein (for plot title)
  l <- dim(data.run2[[1]])[1]
  if(l == 866){p <- "fib.a"}else
    if(l == 491){p <- "fib.b"}else
      if(l == 453){p <- "fib.g"}else
        if(l == 353){p <- "ra33"}else
          if(l == 466){p <- "vim"}else
            if(l == 663){p <- "pad4"}else
            {p <- "unknown protein"}
  
  dfs <- vector(mode = "list", length = length(data))
  for(i in 1:length(data)){
    dat <- data[[i]]
    zero.idx <- which(dat[,1] == 0)
    dat <- dat[-zero.idx,]
    
    df <- try(as.data.frame.matrix(dat), silent = TRUE)
    if(typeof(df) == "character"){ df <- data.frame(dat)}
    dfs[[i]] <- df
  }
  
  if (LEN == "THREE"){
    p0 <- ggplot(dfs[[1]], aes(x = `Total Abd`, y = `Percent Citrullinated`)) + 
      ylab("% Citrullinated") + xlab("log(Abd)") + 
      labs(title = paste0(p, " ", specifier )) + geom_point( color = "black", shape = 1, fill = "white") + 
      geom_point(data = dfs[[2]], color = "red") + geom_point(data = dfs[[3]], color = "black") + 
      theme_classic() + theme(panel.grid.major = element_line(colour = "light grey")) +  scale_x_log10() 
  }else if (LEN == "TWO"){
    p0 <- ggplot(dfs[[1]], aes(x = `Total Abd`, y = `Percent Citrullinated`)) + 
      ylab("% Citrullinated") + xlab("log(Abd)") + 
      labs(title = paste0(p, " ", specifier )) + geom_point( color = "black", shape = 1, fill = "white") + 
      geom_point(data = dfs[[2]], color = "black")  + 
      theme_classic() + theme(panel.grid.major = element_line(colour = "light grey")) +  scale_x_log10() 
  }
  
  
  return(p0)
}
plot.abd.hist <- function(dat, cond, specifier){
  #Input is an r x 2 matrix, with column 1 names "Percent Citrullinated"
  #and column 2 "Total Abd"
  #specifier is a string to add to title. 

  if(missing(cond)){
    cond <- NULL
  }
  
  if(missing(specifier)){
    specifier<- NULL
  }else{
    specifier <- paste0(": ",specifier )
  }
  #Determine protein (for plot title)
  l <- dim(dat)[1]
  if(l == 866){p <- "fib.a"}else
    if(l == 491){p <- "fib.b"}else
      if(l == 453){p <- "fib.g"}else
        if(l == 353){p <- "ra33"}else
          if(l == 466){p <- "vim"}else
            if(l == 663){p <- "pad4"}else
            {p <- "unknown protein"}
  
  df <- as.data.frame.matrix(dat)
  p <- ggplot(df, aes(x = `Total Abd`)) + ylab("# Residues") + xlab("Abundance") + 
    labs(title = paste0(p, " ", cond, specifier )) + geom_histogram(aes()) + theme_minimal() + 
    scale_x_log10()
  return(p)
}