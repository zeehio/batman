createPuerSpectraTemplate <- function(pureSpecDir, metaNames, saveDir)
{
	if (length(pureSpecDir) != length(metaNames))
	{
		stop("Length of PureSpectra and metaNames does not match.\n")
	}
	
	dirPST<-paste(saveDir,"/PureSpectraTemplate",sep="")
	if(!file.exists(dirPST)){
		dir.create(dirPST)
	}
	
	for (n in 1:length(pureSpecDir))
		 {
		 ps <- readBruker(pureSpecDir[n])
		 if (dim(ps)[2] !=2)
		 {
		 stop(paste("More than one spectrum found from row ", n, " of PureSpectra.\n", sep =""))
		 }
		 write.table(ps,file=paste(dirPST,"/",metaNames[n],".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep = "\t")
		 }
}