checkBatmanOptions<-function(dir)
{
  ## written by Dr. Jie Hao, Imperial College London
  dirTime <-NULL
  pBI_IP <- FALSE
  dirsys <- paste(system.file("extdata",package="batman"), '/batmanOptions.txt', sep = '')
  contentsys <-read.table(dirsys,sep="\n",comment.char = "")
  iscommentsys <- substr(contentsys[,1],1,1)
  # comment lines
  cidsys <- which(iscommentsys == '%')
  # parameter lines
  cidnsys <- which(iscommentsys != '%')
  cosys <- as.matrix(contentsys)
  
  
  content <-read.table(dir,sep="\n",comment.char = "")
  iscomment <- substr(content[,1],1,1)
  ## remove empty line
  sID <- NULL
  for (i in length(iscomment):1)
  {
    if (iscomment[i] == ' ')
    {
      sID <- c(sID, i)  
    }      
  }
  iscomment <- iscomment[setdiff(1:length(iscomment), sID)]
  # comment lines
  cid <- which(iscomment == '%')
  # parameter lines
  cidn <- which(iscomment != '%')
  
  content <- content[setdiff(1:dim(content)[1], sID),]
  co <- as.matrix(content)
  
  firsts <- substr(co[cidn,1],1,6)
  rid <- which(firsts == 'rdelta')
  
  cend <- which(cid > cidn[rid])
  
  if (length(cidn) == (length(cidnsys)-4))
  {  
    ## add last line
    co<- rbind(co, 'csFlag - Specify chemical shift for each multiplet in each spectrum? (chemShiftperSpectra.csv file) (Yes - 1 / No - 0): 0')
    pBI_IP <- TRUE
    iscommentco <- substr(co[,1],1,1)
    cidn <- which(iscommentco != '%')
  } 
  if (length(cidn) == (length(cidnsys)-3))
  {
    co2 <- co
    for (i in 1:7)
    {
      co2 <- rbind(co2, co[1,])
    }
    
    j <- 3
    co2[(cidn[j]+1):(dim(co)[1]+1),] <- co[cidn[j]:dim(co)[1],]
    co2[cidn[j],1] <- c('paraProc - No of parallel processes (multicores) if multi-spectra inputted in specNo (only 1 core will be used for single spectrum): 1')  
    
    j <- 9
    co2[(cidn[j]+7):(dim(co)[1]+7),] <- co[cidn[j]:dim(co)[1],]  
    
    
    co2[cidn[j]+1,1] <- c('nItPostBurnin - Number of post-burn-in iterations: 100')   
    co2[cidn[9]+2,1] <- c('multFile - Choose template of multiplets file from options below: 1')
    co2[cidn[9]+3,1] <- c('%% 1, The default template of multiplets in multi_data.csv file,')
    co2[cidn[9]+4,1] <- c('%% 2, The user input template of multiplets in multi_data_user.csv file,')
    co2[cidn[9]+5,1] <- c('%% 3, Both the default and user input template of multiplets files.')
    co2[cidn[9]+6,1] <- c('%%')
    
    co <- co2   
  
    pBI_IP <- TRUE
  } else if (length(cidn) == (length(cidnsys)-1))
  {
    co2 <- co
    co2 <- rbind(co2, co[1,])
    
    j <- 3
    co2[(cidn[j]+1):(dim(co)[1]+1),] <- co[cidn[j]:dim(co)[1],]
    co2[cidn[j],1] <- c('paraProc - No of parallel processes (multicores) if multi-spectra inputted in specNo (only 1 core will be used for single spectrum): 1')  
    co <- co2   

    pBI_IP <- TRUE
  }
  
  iscomment <- substr(co[,1],1,1)
  cid <- which(iscomment == '%')
  cidn <- which(iscomment != '%')
  
  if (pBI_IP)
  {
    if (length(cidnsys) == length(cidn))
    {
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
    
    ctime <- format(Sys.time(), "%d_%b_%H_%M_%S")
    dirTime<-paste(substr(dir,1,nchar(dir)-4), "_", ctime, ".txt",sep="")
    file.rename(from = dir, to = dirTime)
    
    write.table(co, file = dir, sep = "\n", row.names = FALSE,col.names = FALSE,quote=FALSE)
  }
  
  return (list(pBI_IP = pBI_IP, dirTime = dirTime))
}