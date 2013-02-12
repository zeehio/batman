checkBatmanOptions<-function(dir = dir1)
{
## written by Dr. Jie Hao
  content <-read.table(dir,sep="\n",comment.char = "")
  
  iscomment <- substr(content[,1],1,1)
  cid <- which(iscomment == '%')
  cidn <- which(iscomment != '%')
  co <- as.matrix(content)
  
  firsts <- substr(co[cidn,1],1,6)
  rid <- which(firsts == 'rdelta')
  
  cend <- which(cid > cidn[rid])
  
  if (firsts[2] != "Spectr" )
  {
    myVector <- strsplit(co[cidn[2],1], ":")
    co[cidn[2],1]<- c(paste('Spectra range to be included: ', myVector[[1]][2], sep = ""))
  }
  if (length(cidn) > rid)
  {
    
    if (firsts[rid+1] != 'Use sp')
    {
      if (length(cend) == 0)
      {
        co[rid+1,1]<- c('Use specified chemical shift for spectra (chemShiftperSpectra.csv) file (1/0): 0')
      }
      else if (cid[cend[length(cend)]]<dim(content)[1])
      {
        co[cid[cend[length(cend)]]+1,1]<- c('Use specified chemical shift for spectra (chemShiftperSpectra.csv) file (1/0): 0')
      }
      else
      {
        co<- rbind(co, 'Use specified chemical shift for spectra (chemShiftperSpectra.csv) file (1/0): 0')
      }
    } 
    
  } else {
    co<- rbind(co, 'Use specified chemical shift for spectra (chemShiftperSpectra.csv) file (1/0): 0')
  }
  write.table(co, file = dir, sep = "\n", row.names = FALSE,col.names = FALSE,quote=FALSE)
}