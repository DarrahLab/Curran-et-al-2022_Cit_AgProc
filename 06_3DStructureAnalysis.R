path = "C:/Users/13174/OneDrive - Johns Hopkins/Documents/Darrah Lab/Pymol_Ashley/Final source files/"
subdir = "SourceData"
library(readxl)
ra33_native = read_xlsx(paste0(path,subdir, "ra33_native.xlsx"))
ra33_pad4 = read_xlsx(paste0(path,subdir, "ra33_pad4.xlsx"))
ra33_native_p2 = read_xlsx(paste0(path,subdir, "ra33_native_pad2.xlsx"))
ra33_pad2 = read_xlsx(paste0(path,subdir, "ra33_pad2.xlsx"))
fibb_native = read_xlsx(paste0(path,subdir, "fibb_native_pad4.xlsx"))
fibb_pad4 = read_xlsx(paste0(path,subdir, "fibb_pad4.xlsx"))
fibb_native_p2 = read_xlsx(paste0(path,subdir, "fibb_pad2_native.xlsx"))
fibb_pad2 = read_xlsx(paste0(path, "fibb_pad2.xlsx"))
vim_native = read_xlsx(paste0(path,subdir,"vimentin_native_pad4.xlsx"))
vim_pad4 = read_xlsx(paste0(path,subdir, "vim_pad4.xlsx"))
vim_native_p2 = read_xlsx(paste0(path,subdir, "vim_native_p2.xlsx")) 
vim_pad2 = read_xlsx(paste0(path,subdir, "vim_pad2.xlsx"))



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

##################### Make alpha carbon dataframes ############################
###################### Start with RA33 -- PAD4 ################################
###############################################################################

alpha_carbs=data.frame(atom_num=numeric(0), atom_type=character(0), residue_n=character(0), x=numeric(0), y=numeric(0), z=numeric(0), dist_to_nearest_cit=numeric(0))
counter=0
for (i in 1:5088){
  if (ra33_native[[i,3]] == 'CA'){
    counter=counter+1
    alpha_carbs[[counter,1]] = ra33_native[[i,2]]
    alpha_carbs[[counter,2]] = ra33_native[[i,3]]
    alpha_carbs[[counter,3]] = ra33_native[[i,4]]
  }
}
ca_counter=0
for (i in 1:5053){
  if (ra33_pad4[[i,3]] == 'CA'){
    ca_counter=ca_counter+1
    if (alpha_carbs[[ca_counter,3]] == ra33_pad4[[i,4]]){ # As long as it's a matched residue...
      alpha_carbs[[ca_counter,4]] = ra33_pad4[[i,7]]
      alpha_carbs[[ca_counter,5]] = ra33_pad4[[i,8]]
      alpha_carbs[[ca_counter,6]] = ra33_pad4[[i,9]]
    } else {
      alpha_carbs[[ca_counter,3]] ='CIT' # If it's not a matched residue, it's a citrulline.
      alpha_carbs[[ca_counter,4]] = ra33_pad4[[i,7]]
      alpha_carbs[[ca_counter,5]] = ra33_pad4[[i,8]]
      alpha_carbs[[ca_counter,6]] = ra33_pad4[[i,9]]
    }
  }
}

for (i in 1:353){
  print(i)
  dist_nearest_cit = 100000
  for (j in 1:353){
    if (alpha_carbs[[j,3]] == 'CIT'){
      distance = sqrt((alpha_carbs[[i,4]] - alpha_carbs[[j,4]])^2 + (alpha_carbs[[i,5]] - alpha_carbs[[j,5]])^2 + (alpha_carbs[[i,6]] - alpha_carbs[[j,6]])^2) # calculate dist to nearest cit in the pad4 sequence
      if (distance < dist_nearest_cit){
        dist_nearest_cit=distance
      }
    }
  }
  alpha_carbs[[i,7]] = dist_nearest_cit
}

dist_mat = matrix(data=NA, nrow=353, ncol=353)
for (i in 1:353){
  print(i)
  x_1=alpha_carbs[[i,4]]
  y_1=alpha_carbs[[i,5]]
  z_1=alpha_carbs[[i,6]]
  for (j in 1:353){
    x_2=alpha_carbs[[j,4]]
    y_2=alpha_carbs[[j,5]]
    z_2=alpha_carbs[[j,6]]
    dist_mat[[i,j]] = sqrt((x_1-x_2)^2 + (y_1-y_2)^2 + (z_1-z_2)^2)
  }
}


###############################################################################
# Now RA33 PAD2
###############################################################################

alpha_carbs_ra33pad2=data.frame(atom_num=numeric(0), atom_type=character(0), residue_n=character(0), x=numeric(0), y=numeric(0), z=numeric(0), dist_to_nearest_cit=numeric(0))
counter=0
for (i in 1:5088){
  if (ra33_native_p2[[i,3]] == 'CA'){
    counter=counter+1
    alpha_carbs_ra33pad2[[counter,1]] = ra33_native_p2[[i,2]]
    alpha_carbs_ra33pad2[[counter,2]] = ra33_native_p2[[i,3]]
    alpha_carbs_ra33pad2[[counter,3]] = ra33_native_p2[[i,4]]
  }
}
ca_counter=0
for (i in 1:5032){
  if (ra33_pad2[[i,3]] == 'CA'){
    ca_counter=ca_counter+1
    if (alpha_carbs_ra33pad2[[ca_counter,3]] == ra33_pad2[[i,4]]){ # As long as it's a matched residue...
      alpha_carbs_ra33pad2[[ca_counter,4]] = ra33_pad2[[i,7]]
      alpha_carbs_ra33pad2[[ca_counter,5]] = ra33_pad2[[i,8]]
      alpha_carbs_ra33pad2[[ca_counter,6]] = ra33_pad2[[i,9]]
    } else {
      alpha_carbs_ra33pad2[[ca_counter,3]] ='CIT' # If it's not a matched residue, it's a citrulline.
      alpha_carbs_ra33pad2[[ca_counter,4]] = ra33_pad2[[i,7]]
      alpha_carbs_ra33pad2[[ca_counter,5]] = ra33_pad2[[i,8]]
      alpha_carbs_ra33pad2[[ca_counter,6]] = ra33_pad2[[i,9]]
    }
  }
}

for (i in 1:353){
  print(i)
  dist_nearest_cit = 100000
  for (j in 1:353){
    if (alpha_carbs_ra33pad2[[j,3]] == 'CIT'){
      distance = sqrt((alpha_carbs_ra33pad2[[i,4]] - alpha_carbs_ra33pad2[[j,4]])^2 + (alpha_carbs_ra33pad2[[i,5]] - alpha_carbs_ra33pad2[[j,5]])^2 + (alpha_carbs_ra33pad2[[i,6]] - alpha_carbs_ra33pad2[[j,6]])^2) # calculate dist to nearest cit in the pad4 sequence
      if (distance < dist_nearest_cit){
        dist_nearest_cit=distance
      }
    }
  }
  alpha_carbs_ra33pad2[[i,7]] = dist_nearest_cit
}

dist_mat_ra33pad2 = matrix(data=NA, nrow=353, ncol=353)
for (i in 1:353){
  print(i)
  x_1=alpha_carbs_ra33pad2[[i,4]]
  y_1=alpha_carbs_ra33pad2[[i,5]]
  z_1=alpha_carbs_ra33pad2[[i,6]]
  for (j in 1:353){
    x_2=alpha_carbs_ra33pad2[[j,4]]
    y_2=alpha_carbs_ra33pad2[[j,5]]
    z_2=alpha_carbs_ra33pad2[[j,6]]
    dist_mat_ra33pad2[[i,j]] = sqrt((x_1-x_2)^2 + (y_1-y_2)^2 + (z_1-z_2)^2)
  }
}

################################################################################
##### Let's do Fib B 
################################################################################

alpha_carbs_fibb=data.frame(atom_num=numeric(0), atom_type=character(0), residue_n=character(0), x_n=numeric(0), y=numeric(0), z=numeric(0), dist_to_nearest_cit=numeric(0))
counter=0
for (i in 1:7723){
  if (fibb_native[[i,3]] == 'CA'){
    counter=counter+1
    alpha_carbs_fibb[[counter,1]] = fibb_native[[i,2]]
    alpha_carbs_fibb[[counter,2]] = fibb_native[[i,3]]
    alpha_carbs_fibb[[counter,3]] = fibb_native[[i,4]]
  }
}
ca_counter=0
for (i in 1:7695){
  if (fibb_pad4[[i,3]] == 'CA'){
    ca_counter=ca_counter+1
    if (alpha_carbs_fibb[[ca_counter,3]] == fibb_pad4[[i,4]]){ # As long as it's a matched residue
      alpha_carbs_fibb[[ca_counter,4]] = fibb_pad4[[i,7]]
      alpha_carbs_fibb[[ca_counter,5]] = fibb_pad4[[i,8]]
      alpha_carbs_fibb[[ca_counter,6]] = fibb_pad4[[i,9]]
    } else {
      alpha_carbs_fibb[[ca_counter,3]] ='CIT'
      alpha_carbs_fibb[[ca_counter,4]] = fibb_pad4[[i,7]]
      alpha_carbs_fibb[[ca_counter,5]] = fibb_pad4[[i,8]]
      alpha_carbs_fibb[[ca_counter,6]] = fibb_pad4[[i,9]]
    }
  }
}

for (i in 1:491){
  print(i)
  dist_nearest_cit = 100000
  for (j in 1:491){
    if (alpha_carbs_fibb[[j,3]] == 'CIT'){
      distance = sqrt((alpha_carbs_fibb[[i,4]] - alpha_carbs_fibb[[j,4]])^2 + (alpha_carbs_fibb[[i,5]] - alpha_carbs_fibb[[j,5]])^2 + (alpha_carbs_fibb[[i,6]] - alpha_carbs_fibb[[j,6]])^2) 
      if (distance < dist_nearest_cit){
        dist_nearest_cit=distance
      }
    }
  }
  alpha_carbs_fibb[[i,7]] = dist_nearest_cit
}

dist_mat_fibb = matrix(data=NA, nrow=491, ncol=491)

for (i in 1:491){
  print(i)
  x_1=alpha_carbs_fibb[[i,4]]
  y_1=alpha_carbs_fibb[[i,5]]
  z_1=alpha_carbs_fibb[[i,6]]
  for (j in 1:491){
    x_2=alpha_carbs_fibb[[j,4]]
    y_2=alpha_carbs_fibb[[j,5]]
    z_2=alpha_carbs_fibb[[j,6]]
    dist_mat_fibb[[i,j]] = sqrt((x_1-x_2)^2 + (y_1-y_2)^2 + (z_1-z_2)^2)
  }
}


###############################################################################
### Fibb PAD2

alpha_carbs_fibb2=data.frame(atom_num=numeric(0), atom_type=character(0), residue_n=character(0), x = numeric(0), y=numeric(0), z=numeric(0), dist_to_nearest_cit=numeric(0))
counter=0
for (i in 1:7723){
  if (fibb_native_p2[[i,3]] == 'CA'){
    counter=counter+1
    alpha_carbs_fibb2[[counter,1]] = fibb_native_p2[[i,2]]
    alpha_carbs_fibb2[[counter,2]] = fibb_native_p2[[i,3]]
    alpha_carbs_fibb2[[counter,3]] = fibb_native_p2[[i,4]]
  }
}
ca_counter=0
for (i in 1:7688){
  if (fibb_pad2[[i,3]] == 'CA'){
    ca_counter=ca_counter+1
    if (alpha_carbs_fibb2[[ca_counter,3]] == fibb_pad2[[i,4]]){ # As long as it's a matched residue
      alpha_carbs_fibb2[[ca_counter,4]] = fibb_pad2[[i,7]]
      alpha_carbs_fibb2[[ca_counter,5]] = fibb_pad2[[i,8]]
      alpha_carbs_fibb2[[ca_counter,6]] = fibb_pad2[[i,9]]
    } else {
      alpha_carbs_fibb2[[ca_counter,3]] ='CIT'
      alpha_carbs_fibb2[[ca_counter,4]] = fibb_pad2[[i,7]]
      alpha_carbs_fibb2[[ca_counter,5]] = fibb_pad2[[i,8]]
      alpha_carbs_fibb2[[ca_counter,6]] = fibb_pad2[[i,9]]
    }
  }
}
for (i in 1:491){
  print(i)
  dist_nearest_cit = 100000
  for (j in 1:491){
    if (alpha_carbs_fibb2[[j,3]] == 'CIT'){
      distance = sqrt((alpha_carbs_fibb2[[i,4]] - alpha_carbs_fibb2[[j,4]])^2 + (alpha_carbs_fibb2[[i,5]] - alpha_carbs_fibb2[[j,5]])^2 + (alpha_carbs_fibb2[[i,6]] - alpha_carbs_fibb2[[j,6]])^2) # calculate dist to nearest glu in the pad4 sequence
      if (distance < dist_nearest_cit){
        dist_nearest_cit=distance
      }
    }
  }
  alpha_carbs_fibb2[[i,7]] = dist_nearest_cit
}

dist_mat_fibb2 = matrix(data=NA, nrow=491, ncol=491)

for (i in 1:491){
  print(i)
  x_1=alpha_carbs_fibb2[[i,4]]
  y_1=alpha_carbs_fibb2[[i,5]]
  z_1=alpha_carbs_fibb2[[i,6]]
  for (j in 1:491){
    x_2=alpha_carbs_fibb2[[j,4]]
    y_2=alpha_carbs_fibb2[[j,5]]
    z_2=alpha_carbs_fibb2[[j,6]]
    dist_mat_fibb2[[i,j]] = sqrt((x_1-x_2)^2 + (y_1-y_2)^2 + (z_1-z_2)^2)
  }
}

###############################################################################
### Last, we do Vimentin -- PAD4 first
###############################################################################


alpha_carbs_vim=data.frame(atom_num=numeric(0), atom_type=character(0), residue_n=character(0), x=numeric(0), y=numeric(0), z=numeric(0), dist_to_nearest_cit=numeric(0))
counter=0
for (i in 1:7478){
  if (vim_native[[i,3]] == 'CA'){
    counter=counter+1
    alpha_carbs_vim[[counter,1]] = vim_native[[i,2]]
    alpha_carbs_vim[[counter,2]] = vim_native[[i,3]]
    alpha_carbs_vim[[counter,3]] = vim_native[[i,4]]
  }
}
ca_counter=0
for (i in 1:7268){
  if (vim_pad4[[i,3]] == 'CA'){
    ca_counter=ca_counter+1
    if (alpha_carbs_vim[[ca_counter,3]] == vim_pad4[[i,4]]){ # As long as it's a matched residue
      alpha_carbs_vim[[ca_counter,4]] = vim_pad4[[i,7]]
      alpha_carbs_vim[[ca_counter,5]] = vim_pad4[[i,8]]
      alpha_carbs_vim[[ca_counter,6]] = vim_pad4[[i,9]]
    } else {
      alpha_carbs_vim[[ca_counter,3]] ='CIT'
      alpha_carbs_vim[[ca_counter,4]] = vim_pad4[[i,7]]
      alpha_carbs_vim[[ca_counter,5]] = vim_pad4[[i,8]]
      alpha_carbs_vim[[ca_counter,6]] = vim_pad4[[i,9]]
    }
  }
}

for (i in 1:466){
  print(i)
  dist_nearest_cit = 100000
  for (j in 1:466){
    if (alpha_carbs_vim[[j,3]] == 'CIT'){
      distance = sqrt((alpha_carbs_vim[[i,4]] - alpha_carbs_vim[[j,4]])^2 + (alpha_carbs_vim[[i,5]] - alpha_carbs_vim[[j,5]])^2 + (alpha_carbs_vim[[i,6]] - alpha_carbs_vim[[j,6]])^2)
      if (distance < dist_nearest_cit){
        dist_nearest_cit=distance
      }
    }
  }
  alpha_carbs_vim[[i,7]] = dist_nearest_cit
}

dist_mat_vim = matrix(data=NA, nrow=466, ncol=466)
for (i in 1:466){
  print(i)
  x_1=alpha_carbs_vim[[i,4]]
  y_1=alpha_carbs_vim[[i,5]]
  z_1=alpha_carbs_vim[[i,6]]
  for (j in 1:466){
    x_2=alpha_carbs_vim[[j,4]]
    y_2=alpha_carbs_vim[[j,5]]
    z_2=alpha_carbs_vim[[j,6]]
    dist_mat_vim[[i,j]] = sqrt((x_1-x_2)^2 + (y_1-y_2)^2 + (z_1-z_2)^2)
  }
}


################### Finally, Vim PAD2 ##########################################

alpha_carbs_vim_2=data.frame(atom_num=numeric(0), atom_type=character(0), residue_n=character(0), x = numeric(0), y=numeric(0), z=numeric(0), dist_to_nearest_cit=numeric(0))
counter=0
for (i in 1:7478){
  if (vim_native_p2[[i,3]] == 'CA'){
    counter=counter+1
    alpha_carbs_vim_2[[counter,1]] = vim_native_p2[[i,2]]
    alpha_carbs_vim_2[[counter,2]] = vim_native_p2[[i,3]]
    alpha_carbs_vim_2[[counter,3]] = vim_native_p2[[i,4]]
  }
}
ca_counter=0
for (i in 1:7296){
  if (vim_pad2[[i,3]] == 'CA'){
    ca_counter=ca_counter+1
    if (alpha_carbs_vim_2[[ca_counter,3]] == vim_pad2[[i,4]]){ # As long as it's a matched residue
      alpha_carbs_vim_2[[ca_counter,4]] = vim_pad2[[i,7]]
      alpha_carbs_vim_2[[ca_counter,5]] = vim_pad2[[i,8]]
      alpha_carbs_vim_2[[ca_counter,6]] = vim_pad2[[i,9]]
    } else {
      alpha_carbs_vim_2[[ca_counter,3]] ='CIT'
      alpha_carbs_vim_2[[ca_counter,4]] = vim_pad2[[i,7]]
      alpha_carbs_vim_2[[ca_counter,5]] = vim_pad2[[i,8]]
      alpha_carbs_vim_2[[ca_counter,6]] = vim_pad2[[i,9]]
    }
  }
}
for (i in 1:466){
  print(i)
  dist_nearest_cit = 100000
  for (j in 1:466){
    if (alpha_carbs_vim_2[[j,3]] == 'CIT'){
      distance = sqrt((alpha_carbs_vim_2[[i,4]] - alpha_carbs_vim_2[[j,4]])^2 + (alpha_carbs_vim_2[[i,5]] - alpha_carbs_vim_2[[j,5]])^2 + (alpha_carbs_vim_2[[i,6]] - alpha_carbs_vim_2[[j,6]])^2) 
      if (distance < dist_nearest_cit){
        dist_nearest_cit=distance
      }
    }
  }
  alpha_carbs_vim_2[[i,7]] = dist_nearest_cit
}

dist_mat_vim_2 = matrix(data=NA, nrow=466, ncol=466)
for (i in 1:466){
  print(i)
  x_1=alpha_carbs_vim_2[[i,4]]
  y_1=alpha_carbs_vim_2[[i,5]]
  z_1=alpha_carbs_vim_2[[i,6]]
  for (j in 1:466){
    x_2=alpha_carbs_vim_2[[j,4]]
    y_2=alpha_carbs_vim_2[[j,5]]
    z_2=alpha_carbs_vim_2[[j,6]]
    dist_mat_vim_2[[i,j]] = sqrt((x_1-x_2)^2 + (y_1-y_2)^2 + (z_1-z_2)^2)
  }
}


############## Now calculate min actual and exp distances ######################


regions_all=brdf ### Alex FIX THIS
regions_all$dist=0
regions_all$exp=0
regions_all$diff=0
set.seed(42)
for (i in 1:nrow(brdf)){ 
  print(i)
  st=regions_all[[i,1]]
  en=regions_all[[i,2]]
  if (regions_all[[i,5]] =='ra33'){
    if (regions_all[[i,6]] == 'pad4'){
      print('ra33')
      # Import actual distance to nearest citrulline
      regions_all[[i,7]] = min(alpha_carbs$dist_to_nearest_cit[st:en])
      ## OK. But now what is the EXPECTED dist to cit?
      # First, generate 2000 random lists of cits
      cits_lists = rcomb(353,5,2000)
      stoch_dists_to_nearest_cit = c()
      for (l in 1:2000){
        cits=cits_lists[l]
        minimum=10000
        for (o in st:en){
          for (p in 1:5){
            row=as.numeric(lapply(cits_lists[l], "[[", p))
            distance=dist_mat[[o,row]]
            if (distance < minimum){
              minimum=distance
            }
          }
        }
        stoch_dists_to_nearest_cit = c(stoch_dists_to_nearest_cit, minimum)
      }
      regions_all[[i,8]] = mean(stoch_dists_to_nearest_cit)
      ################## Now let's do pad2 ##########################
    } else if (regions_all[[i,6]] == "pad2"){
      print('ra33')
      regions_all[[i,7]] = min(alpha_carbs_ra33pad2$dist_to_nearest_cit[st:en])
      cits_lists = rcomb(353,8,2000)
      stoch_dists_to_nearest_cit = c()
      for (l in 1:2000){
        cits=cits_lists[l]
        minimum=10000
        for (o in st:en){
          for (p in 1:8){ # 8 citrullines in pad2-cit ra33
            row=as.numeric(lapply(cits_lists[l], "[[", p))
            distance=dist_mat_ra33pad2[[o,row]]
            if (distance < minimum){
              minimum=distance
            }
          }
        }
        stoch_dists_to_nearest_cit = c(stoch_dists_to_nearest_cit, minimum)
      }
      regions_all[[i,8]] = mean(stoch_dists_to_nearest_cit)
    } 
    # Things that change with a new antigen --
    # 1 - the alpha carbon dataframe
    # 2 - the number of citrullines to generate in the bootstrap
    # 3 - the distance matrix used
    # 4 - the number of alpha carbons to look through
  } else if (regions_all[[i,5]] == "fibb"){
    if (regions_all[[i,6]] == "pad4"){
      print('fibb')
      regions_all[[i,7]] = min(alpha_carbs_fibb$dist_to_nearest_cit[st:en])
      cits_lists = rcomb(353,4,2000) # Only 4 citrullines in fibb pad 4
      stoch_dists_to_nearest_cit = c()
      for (l in 1:2000){
        cits=cits_lists[l]
        minimum=10000
        for (o in st:en){
          for (p in 1:4){
            row=as.numeric(lapply(cits_lists[l], "[[", p))
            distance=dist_mat_fibb[[o,row]]
            if (distance < minimum){
              minimum=distance
            }
          }
        }
        stoch_dists_to_nearest_cit = c(stoch_dists_to_nearest_cit, minimum)
      }
      regions_all[[i,8]] = mean(stoch_dists_to_nearest_cit)
      # Now do the same for pad2
    } else if (regions_all[[i,6]] == "pad2"){
      print('fibb')
      regions_all[[i,7]] = min(alpha_carbs_fibb2$dist_to_nearest_cit[st:en])
      cits_lists = rcomb(353,5,2000) # 4 cits in pad 2 fibb
      stoch_dists_to_nearest_cit = c()
      for (l in 1:2000){
        cits=cits_lists[l]
        minimum=10000
        for (o in st:en){
          for (p in 1:5){
            row=as.numeric(lapply(cits_lists[l], "[[", p))
            distance=dist_mat_fibb2[[o,row]]
            if (distance < minimum){
              minimum=distance
            }
          }
        }
        stoch_dists_to_nearest_cit = c(stoch_dists_to_nearest_cit, minimum)
      }
      regions_all[[i,8]] = mean(stoch_dists_to_nearest_cit)
    }
    
    
  } else if (regions_all[[i,5]] == "vim"){
    if (regions_all[[i,6]] == "pad4"){
      print('vim')
      regions_all[[i,7]] = min(alpha_carbs_vim$dist_to_nearest_cit[st:en])
      cits_lists = rcomb(353,30,2000)
      stoch_dists_to_nearest_cit = c()
      for (l in 1:2000){
        cits=cits_lists[l]
        minimum=10000
        for (o in st:en){
          for (p in 1:30){
            row=as.numeric(lapply(cits_lists[l], "[[", p))
            distance=dist_mat_vim[[o,row]]
            if (distance < minimum){
              minimum=distance
            }
          }
        }
        stoch_dists_to_nearest_cit = c(stoch_dists_to_nearest_cit, minimum)
      }
      regions_all[[i,8]] = mean(stoch_dists_to_nearest_cit)
    } else if (regions_all[[i,6]] == "pad2"){
      print('vim')
      regions_all[[i,7]] = min(alpha_carbs_vim_2$dist_to_nearest_cit[st:en])
      cits_lists = rcomb(353,26,2000) # 26 cits in vim pad2
      stoch_dists_to_nearest_cit = c()
      for (l in 1:2000){
        cits=cits_lists[l]
        minimum=10000
        for (o in st:en){
          for (p in 1:26){
            row=as.numeric(lapply(cits_lists[l], "[[", p))
            distance=dist_mat_vim_2[[o,row]]
            if (distance < minimum){
              minimum=distance
            }
          }
        }
        stoch_dists_to_nearest_cit = c(stoch_dists_to_nearest_cit, minimum)
      }
      regions_all[[i,8]] = mean(stoch_dists_to_nearest_cit)
    }
  }
}

regions_all$diff = regions_all$dist - regions_all$exp
saved_file = save(regions_all, file="regions_all_4_13.rda")


