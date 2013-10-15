checkBatmanOptions<-function(dir)
{
  ## written by Dr. Jie Hao, Imperial College London
  
  pBI_IP <- FALSE
  content <-read.table(dir,sep="\n",comment.char = "")
  
  iscomment <- substr(content[,1],1,1)
  # comment lines
  cid <- which(iscomment == '%')
  # parameter lines
  cidn <- which(iscomment != '%')
  co <- as.matrix(content)
  
  firsts <- substr(co[cidn,1],1,6)
  rid <- which(firsts == 'rdelta')
  
  cend <- which(cid > cidn[rid])
  
  #if (firsts[2] != 'Spectr' )
  #{
    if (firsts[2] != 'specNo')
    {
      myVector <- strsplit(co[cidn[2],1], ":")
      co[cidn[2],1]<- c(paste('specNo - Ranges of spectra number to be included (e.g. 1,3-4 etc.): ', myVector[[1]][2], sep = "")) 
    }
  #}
  if (length(cidn) > rid)
  {    
    #if (firsts[rid+1] != 'Use sp')
    #{
      if (firsts[rid+1] != 'csFlag')
      {
        if (length(cend) == 0)
        {
          co[cidn[rid]+1,1]<- c('csFlag - Specify chemical shift for each multiplet in each spectrum? (chemShiftperSpectra.csv file) (Yes - 1 / No - 0): 0')
        }
        else if (cid[cend[length(cend)]]<dim(content)[1])
        {
          co[cid[cend[length(cend)]]+1,1]<- c('csFlag - Specify chemical shift for each multiplet in each spectrum? (chemShiftperSpectra.csv file) (Yes - 1 / No - 0): 0')
        }
        else
        {
          co<- rbind(co, 'csFlag - Specify chemical shift for each multiplet in each spectrum? (chemShiftperSpectra.csv file) (Yes - 1 / No - 0): 0')
        }
      } 
    #}
  } else {
    co<- rbind(co, 'csFlag - Specify chemical shift for each multiplet in each spectrum? (chemShiftperSpectra.csv file) (Yes - 1 / No - 0): 0')
  }
  
  #nItPostBurnin - Number of post-burn-in iterations: 100
  if (firsts[9] != 'nItPos' && firsts[8] == 'nItBur')
  {
    co2 <- co
    for (i in 1:6)
    {
      co2 <- rbind(co2, co[1,])
    }
    
    co2[cidn[9],1] <- c('nItPostBurnin - Number of post-burn-in iterations: 100')    
    co2[cidn[9]+1,1] <- c('multFile - Choose template of multiplets file from options below: 2')
    co2[cidn[9]+2,1] <- c('%% 1, The default template of multiplets in multi_data.csv file,')
    co2[cidn[9]+3,1] <- c('%% 2, The user input template of multiplets in multi_data_user.csv file,')
    co2[cidn[9]+4,1] <- c('%% 3, Both the default and user input template of multiplets files.')
    co2[cidn[9]+5,1] <- c('%%')
    co2[(cidn[9]+6):dim(co2)[1],1] <- co[cidn[9]:dim(co)[1],1] 
    co <- co2    
    pBI_IP <- TRUE
  }
  
  
  if (firsts[1] != 'ppmRan')
  {
    dirsys <- paste(system.file("extdata",package="batman"), '/batmanOptions.txt', sep = '')
    contentsys <-read.table(dirsys,sep="\n",comment.char = "")
    
    iscommentsys <- substr(contentsys[,1],1,1)
    # comment lines
    cidsys <- which(iscommentsys == '%')
    # parameter lines
    cidnsys <- which(iscommentsys != '%')
    cosys <- as.matrix(contentsys)
    
    if (length(cidnsys) == length(cidn))
    {
    for (i in 1:length(cidnsys))
    {
      myV1 <- strsplit(cosys[cidnsys[i],1], ":")
      myV2 <- strsplit(co[cidn[i],1], ":")
      cosys[cidnsys[i],1] <- c(paste(myV1[[1]][1], ':', myV2[[1]][2], sep = ""))
    }  
    }
    if (cid[length(cid)] > cidn[length(cidn)])
    {
      for (j in (cidn[length(cidnsys)]+1):length(co))
      cosys<- rbind(cosys, co[j])      
    }
    co <- cosys
  }
  write.table(co, file = dir, sep = "\n", row.names = FALSE,col.names = FALSE,quote=FALSE)
  return (pBI_IP)
}