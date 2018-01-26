saveBruker2Txt <- function(BrukerDataDir,saveFileName, sep = "\t")
{
  ## written by Dr. Jie Hao, Imperial College London
  sa<-readBrukerZipped(BrukerDataDir)  
  write.table(sa,file=saveFileName,row.names=FALSE,col.names=TRUE,quote=FALSE,sep = sep)
}
