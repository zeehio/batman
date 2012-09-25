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
  dirTPure<-paste(dirIP,"/testBrukerPureSpec",sep="")
  if(!file.exists(dirTPure)){
		dir.create(dirTPure)
	}
  
  dirTPure1<-paste(dirTPure,"/testPure1",sep="")
  if(!file.exists(dirTPure1)){
    dir.create(dirTPure1)
  }
  
  dirTPure2<-paste(dirTPure,"/testPure2",sep="")
  if(!file.exists(dirTPure2)){
    dir.create(dirTPure2)
  }
  
  dirOP<-paste(dirRoot,"/BatmanOutput",sep="")
  if(!file.exists(dirOP)) {
    dir.create(dirOP)
  }
  extdir<-system.file("extdata",package="batman")
	
	src.dir <- paste(extdir,"/testBrukerPureSpec/testPure1",sep="")
	dest.dir <-dirTPure1
	
	file.names <- dir(src.dir)
	sapply(file.names, function(x) {
		   file.copy(from=paste(src.dir, x, sep="/"),
					 to=paste(dest.dir, x, sep="/"),
					 overwrite = FALSE) })
	
	src.dir <- paste(extdir,"/testBrukerPureSpec/testPure2",sep="")
	dest.dir <-dirTPure2
	
	file.names <- dir(src.dir)
	sapply(file.names, function(x) {
		   file.copy(from=paste(src.dir, x, sep="/"),
					 to=paste(dest.dir, x, sep="/"),
					 overwrite = FALSE) })
	
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
  
  dirA<-c(dirRoot, dirIP, dirOP, dirPST, dirTPure)
  ##, dirPST
  return (dirA)
}