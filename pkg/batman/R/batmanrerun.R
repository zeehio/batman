batmanrerun<- function(BM, figBatmanFit = TRUE, listMeta = FALSE, 
figRelCon = FALSE, figMetaFit = FALSE)
{   
## written by Dr. Jie Hao, Imperial College London
## rerun batman() with fixed ppm position for all multiplets
## directories for input parameters and data files
	dir1<-paste(BM$outputDir,"/batmanOptions.txt",sep="")
	dir2<-paste(BM$outputDir,"/NMRdata.txt",sep="")
	dir3<-paste(BM$outputDir,"/metabolitesList.txt",sep="")
	dir4<-paste(BM$outputDir,"/multi_data.dat",sep="")
	dir6<-paste(BM$outputDir,"/chemShiftPerSpec.dat",sep="")
	dir7<-paste(BM$inputDir,"/PureSpectraTemplate/",sep = "")

	warnDef<-options("warn")$warn
	warnRead<-options(warn = -1)
	
## check for duplicated list metabolites
	dirL<-paste(BM$outputDir,"/metabolitesList.csv",sep="")
	mL<-read.csv(dirL, header=F,colClasses="character")
	if (anyDuplicated(mL[,1])!=0)
    stop("Duplicated metabolite list.\n")
	
	mL<-mL[,1,drop=FALSE]
	write.table(mL,file=dir3, sep = "\t", row.names=FALSE,col.names=FALSE,quote=FALSE)
	
	dir5<-paste(BM$outputDir,"/",sep="")
## file dirs for input parameters and data.
	filedir<-c(dir1,dir2,dir3,dir4,dir5,dir6,dir7)
	
## read in parameters from batmanOption.txt 
	nlines <-dim(read.table(dir1,sep="\n",comment.char = ""))[1]
	con  <- file(dir1, open = "r")
	oneLine <- readLines(con, n = nlines, warn = FALSE)
	fL<-substr(oneLine,1,1)
	nL<-which(is.na(match(fL,"%")))
	myVector <- strsplit(oneLine[nL[2]], ":")
	sno <- NULL
	sno <-getSpectraRange(myVector)
	myVector <- strsplit(oneLine[nL[24+1]], ":")
	chemshif <- as.numeric(myVector[[1]][2])
	myVector <- strsplit(oneLine[nL[8+1]], ":")
	itoBI <- as.numeric(myVector[[1]][2])
	myVector <- strsplit(oneLine[nL[10+2+1]], ":")
	fixeff <- as.numeric(myVector[[1]][2])
	myVector <- strsplit(oneLine[nL[11+2+1]], ":")
	itoRr <- as.numeric(myVector[[1]][2])
	## parallel process
	myVector <- strsplit(oneLine[nL[3]], ":") 
	paraProc <- as.numeric(myVector[[1]][2])
	close(con)
	
	if (chemshif == 1)
	{
		if (!file.exists(dir6))
		{
## createChemShift(tempOption = opt, tempData = b, sTitle = stit, dirA[2])
			stop("No chemshifPerSpec.dat file found in BatmanOutput folder.\n
				 Please check batmanOptions.txt for chemShift setting.\n")
		} 
	}
	
	
	cat("Rerunning batman for ", itoRr, " iterations.\n")
## choose whether to parallelize between spectra
	if (length(sno)>1 && fixeff == 0)    
	{
	  wr<- paraProc 
	  cat("\nNumber of parallel processes (multicores) used to run the multi-spectra analysis: ", wr, "\n")
	  cat("\n Parallel processing of multi spectra currently cannot display the progress\n")
	  cat(" bar (or any words), please be patient for the results :)\n\n")    
	} else {
		wr<-1
	}
	b<-read.delim(dir4)
	bn<-nrow(b)
	
## running MCMC in rerun iterations,
	rr <- 1
	if (rr == 1) {       
		if (wr>1) {
			cl<-makeCluster(wr, type = "SOCK")
			registerDoSNOW(cl) 
		} else {
			cat ("percentage completed...\n")
		}
		if (wr>1) {
			stime1 <- system.time ({
								   multispec1<- foreach (nospec = sno,.packages="batman") %dopar% {
								   pBar <- txtProgressBar(min =0, max = itoRr, style = 3)
								   out<-.Call("batman",filedir,as.integer(bn),as.integer(itoRr),as.integer(rr),as.integer(nospec-1),pBar,PACKAGE = "batman")
								   close( pBar )
								   }})[3]} else if (wr == 1 && fixeff == 1) {
				stime1 <- system.time ({
									   multispec1<- for(nospec in 1:1) {
									   pBar <- txtProgressBar(min =0, max = itoRr, style = 3)
									   out<-.Call("batman",filedir,as.integer(bn),as.integer(itoRr),as.integer(rr),as.integer(length(sno)-1),pBar,PACKAGE = "batman")
									   close( pBar ) 
									   }})[3]}  else {
					stime1 <- system.time ({
										   multispec1<- for(nospec in sno) {
										   pBar <- txtProgressBar(min =0, max = itoRr, style = 3)
										   out<-.Call("batman",filedir,as.integer(bn),as.integer(itoRr),as.integer(rr),as.integer(nospec-1),pBar,PACKAGE = "batman")
										   close( pBar ) 
										   }})[3]}	
		cat (" For batman rerun, time ")
		print(stime1)
		cat ("seconds.\n")
		if (wr>1) {stopCluster(cl)}
	}
	
	if (out[2] == 1)
    return(cat("No metabolites with resonances in the region for analysis, exiting ...\n"))
	
	cat("Reading in saved data in folder\n")
## read in results into R
	BMR<-readBatmanOutput(BM$outputDir, BM$inputDir) 
	cat(BMR$outputDir)
## plotting results
	if (figBatmanFit)
    plotBatmanFit(BMR, saveFigDir = BM$outputDir, listMeta = listMeta, rerun = TRUE)
	if (figRelCon) 
    plotRelCon(BMR, saveFigDir = BM$outputDir, rerun = TRUE)
	if (figMetaFit)
    plotMetaFit(BMR, saveFigDir = BM$outputDir, rerun = TRUE)
	
	file.remove(dir3)
	cat("\nCompleted.\n")
	warnRead<-options(warn = warnDef)
	return(BMR)
}
