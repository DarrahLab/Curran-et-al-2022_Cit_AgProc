#Cit-Abd Corr functions. 

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
citNeigbors <- function(abd, d,cit_sites){
  #Abd gives abundances by residue. 
  #d is the distance +/- a residue considered local. 
  #cit_sites is a vector indicating cit sites. 
  
  
  nCit <- length(cit_sites)
  #Determine Neighborhoods 
  #Store neighborhood abd as rows in matrix, ignoring out-of-bounds
  #compute av
  #output result 
}
d_to_cit <- function(res,cit_sites,dir_arg,AgLen){
  #For each residue of interest specified in res, determine distance to cit. 
  #res marks the residues of interest. 
  #dir_arg: 0 or NULL if bidirectional shortest dist, -1 if only searching behind, 
  #         +1 if only searching forward in sequence. 
  
  #Test case: 
  # cit_sites <- c(8,14); res <- c(5,11,14,16); dir_arg <- 0; AgLen <- 20
  # cit_sites <- c(8,14); res <- c(5,11,14,16); dir_arg <- 1; AgLen <- 20
  # cit_sites <- c(8,14); res <- c(5,11,14,16); dir_arg <- -1; AgLen <- 20
  
  #Order data if not already 
  cit_sites <- cit_sites[order(cit_sites)] 
  res <- res[order(res)]
  
  dists <- vector(mode = "integer", length = length(res)); dists[] = NA
  
  for(r in 1:length(res)){
    if(dir_arg == -1){
      cits <- cit_sites[which(cit_sites <= res[r])]
      dmin <- res[r] - cit_sites[length(cits)]
    }else if(dir_arg==1){
      cits <- cit_sites[which(cit_sites >= res[r])]
      dmin <- cits[1] - res[r]
    }else{
      dmin <- min(abs(res[r] - cit_sites))
    }
    
    if(length(dmin) > 0){
    dists[r] <- dmin
    } else {dists[r] <- NA}
  }
  output <- list()
  output$mean <- mean(dists, na.rm = TRUE)
  names(dists) <- paste0("r",res)
  output$dists <- dists
  
  return(output)
}
dCit_compare <- function(args, abd, cit_sites, dir_arg){
  #Determine how the dCit of your input args compares to expected. 
  #Expected dCit determined by monte carlo simulation. 
  #args has N rows corresponding to N args, and 2 cols corresponding to start and stop
  #         position of that fragment. 
  #Abd is an abundance vector by residue. 
  #cit_sites is a vector indicating cit sites. 
  #dir_arg: 0 or NULL if bidirectional shortest dist, -1 if only searching behind, 
  #         +1 if only searching forward in sequence. 
  
  nArgs <- dim(args)[1] 
  
  for(f in 1:nArgs){
    
  }
}
process1 <- function(p2,n){
  #Take log2 ratio between raw abundance vectors. 
  #As shown, this is between pad2 and native w/ pad2 enriched residues being positive
  z.n <- which(n == 0)
  z.p2 <- which(p2 == 0)
  p2.n <- p2 / n
  both.z <- z.p2[which(z.p2 %in% z.n)]
  just.n <- z.n[which((z.n %in% z.p2) == FALSE)]
  just.p2 <- z.p2[which((z.p2 %in% z.n) == FALSE)]
  
  p2.n[both.z] <- 1 # 0 / 0 == 1
  if(length(just.n) > 0){p2.n[just.n] = Inf}
  p2.n <- log(p2.n, base = 2)
  if(length(just.p2) > 0){p2.n[just.p2] = -Inf}
  return(p2.n)
}
#Plot0: good ol' abundance
#Plot1: DiffEx between two conditions - centers around zero, two cutoff lines
#Plot2: Same as Plot1 but labels provided cit sites. 
#Plot3: Same as Plot0 but labels provided cit sites. 
#Plot4X: Me bullshitting some patterns for AC 6 Apr 2022. Mods on Plot2 Fn. 
plot0 <- function(d, cutoff, specifier){
  #log2 the data
  d[which(d == 0)] <- 1
  d <- log(d)
  
  if(missing(cutoff)){
    cutoff <- 6
  }
  if(missing(specifier)){
    specifier <- NULL
  }
  t <- paste0(specifier, " Abundances")
  lim <- max(abs(d[abs(d) < Inf]))
  
  
    p <- ggplot(data = data.frame(x = c(1:length(d)), y = d), aes(x,y)) + 
      geom_area() + ylim(0,lim) + 
      labs(title = t, x = "AA Position", y = "Log(Abd)") + theme_minimal() + theme(text = element_text(size = 15)) + 
      geom_hline(yintercept = 0) +
      geom_hline(yintercept = cutoff, color = "red", linetype = "dashed")
  
}
plot1 <- function(d, cutoff, specifier){
  if(missing(cutoff)){
    cutoff <- NULL
  }
  if(missing(specifier)){
    specifier <- NULL
  }
  t <- paste0("Diff Ex: ", specifier)
  lim <- max(abs(d[abs(d) < Inf]))
  
  if(is.null(cutoff)){
  p <- ggplot(data = data.frame(x = c(1:length(d)), y = d), aes(x,y)) + 
       geom_area() + ylim(-lim,lim) + 
       labs(title = t, y = "log2(Abd)") + theme_minimal() + theme(text = element_text(size = 15)) + 
       geom_hline(yintercept = 0) 
  }else{
    p <- ggplot(data = data.frame(x = c(1:length(d)), y = d), aes(x,y)) + 
      geom_area() + ylim(-lim,lim) + 
      labs(title = t) + theme_minimal() + theme(text = element_text(size = 15)) + 
      geom_hline(yintercept = 0) +
      geom_hline(yintercept = cutoff, color = "red", linetype = "dashed") + 
      geom_hline(yintercept = -cutoff, color = "red", linetype = "dashed") 
  }
  
}
plot2 <- function(d, cit_sites, cutoff, specifier){
  if(missing(cutoff)){
    cutoff <- NULL
  }
  if(missing(specifier)){
    specifier <- NULL
  }
  t <- paste0("Diff Ex: ", specifier)
  lim <- max(abs(d[abs(d) < Inf]))
  
  #Make cit_sites df 
  x <- cit_sites
  y <- lim
  cits <- data.frame(x = x, y = y)
  cits2 <- data.frame(x = x, y = -y)
  
  if(is.null(cutoff)){
    p <- ggplot(data = data.frame(x = c(1:length(d)), y = d), aes(x,y)) + 
     ylim(-lim,lim) + 
      geom_point(data = cits, col = "grey", fill = "red", shape = "triangle down filled",  size = 3) + 
      labs(title = t, y = "log2(Abd)") + theme_minimal() + theme(text = element_text(size = 15)) + 
      geom_hline(yintercept = 0, color = "black") 
  }else{
   
    p <- ggplot(data = data.frame(x = c(1:length(d)), y = d), aes(x,y)) + 
      ylim(-lim,lim) +  
      geom_point(data = cits, col = "grey", fill = "red", shape = "triangle down filled",  size = 3) + 
      geom_area() + 
      labs(title = t, y = "log2(Abd)") + theme_minimal() + theme(text = element_text(size = 15)) + 
      geom_hline(yintercept = 0, color = "black") + 
      geom_hline(yintercept = cutoff, color = "red", linetype = "dashed") + 
      geom_hline(yintercept = -cutoff, color = "red", linetype = "dashed") 
  }
  
}
plot3 <- function(d, cit_sites, cutoff, specifier){
  #log2 the data
  d[which(d == 0)] <- 1
  d <- log(d)
  
  if(missing(cutoff)){
    cutoff <- 6
  }
  if(missing(specifier)){
    specifier <- NULL
  }
  t <- paste0(specifier, " Abundances")
  lim <- max(abs(d[abs(d) < Inf]))
  
  #Make cit_sites df 
  x <- cit_sites
  y <- lim
  cits <- data.frame(x = x, y = y)
  
  p <- ggplot(data = data.frame(x = c(1:length(d)), y = d), aes(x,y)) + 
    geom_area() + ylim(0,lim) + 
    labs(title = t, x = "AA Position", y = "Log(Abd)") + theme_minimal() + theme(text = element_text(size = 15)) + 
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = cutoff, color = "red", linetype = "dashed") + 
    geom_point(data = cits, col = "grey", fill = "red", shape = "triangle down filled",  size = 3)
  
}
plot41 <- function(d, cit_sites, cutoff, specifier){
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
      theme_prism() + 
      geom_point(data = cits, col = "black", fill = "white", shape =21,  size = 3) + 
      ylim(-lim,lim)
    
    return(p)

  
}