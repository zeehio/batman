checkBatmanOptions<-function(dir)
{
  ## written by Dr. Jie Hao, Imperial College London
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
  
  if (firsts[2] != 'Spectr' )
  {
    if (firsts[2] != 'specNo')
    {
      myVector <- strsplit(co[cidn[2],1], ":")
      co[cidn[2],1]<- c(paste('specNo - Ranges of spectra number to be included (e.g. 1,3-4 etc.): ', myVector[[1]][2], sep = "")) 
    }
  }
  if (length(cidn) > rid)
  {    
    if (firsts[rid+1] != 'Use sp')
    {
      if (firsts[rid+1] != 'csFlag')
      {
        if (length(cend) == 0)
        {
          co[rid+1,1]<- c('csFlag - Specify chemical shift for each multiplet in each spectrum? (chemShiftperSpectra.csv file) (Yes - 1 / No - 0): 0')
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
    }
  } else {
    co<- rbind(co, 'csFlag - Specify chemical shift for each multiplet in each spectrum? (chemShiftperSpectra.csv file) (Yes - 1 / No - 0): 0')
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
    
    for (i in 1:length(cidnsys))
    {
      myV1 <- strsplit(cosys[cidnsys[i],1], ":")
      myV2 <- strsplit(co[cidn[i],1], ":")
      cosys[cidnsys[i],1] <- c(paste(myV1[[1]][1], ':', myV2[[1]][2], sep = ""))
    }  
    if (cid[length(cid)] > cidn[length(cidn)])
    {
      for (j in (cidn[length(cidnsys)]+1):length(co))
      cosys<- rbind(cosys, co[j])      
    }
    co <- cosys
  }
  write.table(co, file = dir, sep = "\n", row.names = FALSE,col.names = FALSE,quote=FALSE)
}