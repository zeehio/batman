createPureSpectraTemplate <- function(dirPureSpec, metaNames, dirIP)
{
  if (length(dirPureSpec) != length(metaNames))
  {
    stop("Length of PureSpectra and metaNames does not match.\n")
  }
  
  dirPST<-paste(dirIP,"/PureSpectraTemplate",sep="")
  if(!file.exists(dirPST)){
    dir.create(dirPST)
  }
  
  for (n in 1:length(dirPureSpec))
  {
    ps <- readBruker(dirPureSpec[n])
    if (dim(ps)[2] !=2)
    {
      stop(paste("More than one spectrum found from row ", n, " of PureSpectra:\n", dirPureSpec[n], "\n", sep =""))
    }
    if (ps[1,1] > ps[2,1])
    {
      ps <- ps[nrow(ps):1,]
    }
    write.table(ps,file=paste(dirPST,"/",metaNames[n],".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep = "\t")
  }
}