newDir<-function(runBATMANDir = getwd(), overwriteFile = FALSE)
{
  ## create a new directory with installed data files
  dirRoot<-paste(runBATMANDir,"/runBATMAN",sep="")
  if(!file.exists(dirRoot)) {
    dir.create(dirRoot)
  }
  dirIP<-paste(dirRoot,"/BatmanInput",sep="")
  if(!file.exists(dirIP)){
    dir.create(dirIP)
  }
  dirPST<-paste(dirIP,"/PureSpectraTemplate",sep="")
  if(!file.exists(dirPST)){
    dir.create(dirPST)
  }
  dirOP<-paste(dirRoot,"/BatmanOutput",sep="")
  if(!file.exists(dirOP)) {
    dir.create(dirOP)
  }
  extdir<-system.file("extdata",package="batman")
  dirbo<-paste(extdir,"/batmanOptions.txt",sep="")
  cp<-file.copy(dirbo,dirIP, overwrite=overwriteFile)
  dird<-paste(extdir,"/NMRdata.txt",sep="")
  cp<-file.copy(dird,dirIP, overwrite=overwriteFile)
  dirL<-paste(extdir,"/metabolitesList.csv",sep="")
  cp<-file.copy(dirL,dirIP, overwrite=overwriteFile)
  dirR<-paste(extdir,"/multi_data.csv",sep="")
  cp<-file.copy(dirR,dirIP,overwrite=overwriteFile)
  dirRU<-paste(extdir,"/multi_data_user.csv",sep="")
  cp<-file.copy(dirRU,dirIP,overwrite=overwriteFile)
  
  dirA<-c(dirRoot, dirIP, dirOP,dirPST)
  ##, dirPST
  return (dirA)
}