batman<-function(BrukerDataDir, txtFile, rData, createDir = TRUE, runBATMANDir = getwd(), 
                 overwriteDir = FALSE, figBatmanFit = TRUE, listMeta = FALSE, 
                 figRelCon = FALSE, figMetaFit = FALSE)
{
  warnDef<-options("warn")$warn
  warnRead<-options(warn = -1)
  ## main function
  ## input data files directory 
  if (createDir) {
    dirA<-newDir(runBATMANDir = runBATMANDir, overwriteFile = overwriteDir)
  } else {
    extdir<-system.file("extdata",package="batman")
    dirPST<-paste(extdir,"/PureSpectraTemplate",sep="")
    if(!file.exists(dirPST)){
      dir.create(dirPST)
      dirA<-c(extdir,extdir,extdir, dirPST)
    }
  }
  ctime <- format(Sys.time(), "%d_%b_%H_%M_%S")
  
  dir1<-paste(dirA[2],"/batmanOptions.txt",sep="")
  dir2<-paste(dirA[2],"/NMRdata.txt",sep="")
  dir3<-paste(dirA[2],"/metabolitesList.txt",sep="")
  dir4<-paste(dirA[2],"/multi_data.dat",sep="")
  dir6<-paste(dirA[2],"/chemShiftPerSpec.dat",sep="")
  dir7<-paste(dirA[4],"/",sep = "")
  
  dirctime<-paste(dirA[3],"/",ctime,sep="")
  if(!file.exists(dirctime)) {
    dir.create(dirctime)
  }
  
  dir5<-paste(dirctime,"/",sep="")
  
  cp<-file.copy(dir1,dir5)
  dirList<-paste(dirA[2],"/metabolitesList.csv",sep="")
  cp<-file.copy(dirList,dir5)
  
  filedir<-c(dir1,dir2,dir3,dir4,dir5,dir6,dir7)
  
  dirL<-paste(dirA[2],"/metabolitesList.csv",sep="")
  dirR<-paste(dirA[2],"/multi_data.csv",sep="")
  dirRU<-paste(dirA[2],"/multi_data_user.csv",sep="")
  
  ## read in batman optitons
  con  <- file(dir1, open = "r")
  oneLine <- readLines(con, n = 49, warn = FALSE)
  fL<-substr(oneLine,1,1)
  nL<-which(is.na(match(fL,"%")))
  myVector <- strsplit(oneLine[nL[2]], ":")
  sno <- NULL
  sno <-getSpectraRange(myVector)
  
  myVector <- strsplit(oneLine[nL[24]], ":")
  chemshif <- as.numeric(myVector[[1]][2])
  myVector <- strsplit(oneLine[nL[8]], ":")
  itoBI <- as.numeric(myVector[[1]][2])
  myVector <- strsplit(oneLine[nL[10]], ":")
  fixeff <- as.numeric(myVector[[1]][2])
  myVector <- strsplit(oneLine[nL[11]], ":")
  itoRr <- as.numeric(myVector[[1]][2])
  close(con)
  
  cat("batman...\n")
  cat("Enter number of post-burn-in iterations (burn-in currently set to ",itoBI, " iterations):\n" )
  itoPBI<- getinput(lowlim=1,highlim=-1)  
  
  ito<-itoPBI+itoBI
  ## choose template file
  choices<-c("Include the default template of multiplets in multi_data.csv file only.",
             "Include the user input template of multiplets in multi_data_user.csv file only.",
             "Include both the above files.")
  showLine <-"\nEnter a number of choice from the menu below:\n"
  opt <- menuA(choices, 1, showLine)
  
  if (opt == 1)
  {
    cat("Loading multi_data.csv...\n")
    b<-read.csv(dirR,header=T,stringsAsFactors=FALSE,colClasses="character")
  } else if (opt == 2)  {
    cat("Loading multi_data_user.csv...\n")
    b<-read.csv(dirRU,header=T,stringsAsFactors=FALSE,colClasses="character")
  } else if (opt == 3)  {
    cat("Loading multi_data.csv and multi_data_user.csv...\n")
    b1<-read.csv(dirR,header=T,stringsAsFactors=FALSE,colClasses="character")
    b2<-read.csv(dirRU,header=T,stringsAsFactors=FALSE,colClasses="character")
    b<-rbind(b1,b2)
  } 
  ## prepare template file for c++
  bn<-nrow(b)
  ind <-6:8
  for (i in ind)
  {
    check<-b[,i]=='n'
    b[check,i]=-50
  }    
  D<-order(b[,1])
  bor<-b[D,]
  write.table(bor,file=dir4,sep = "\t",row.names = FALSE,col.names = FALSE,quote=FALSE)
  
  mL<-read.csv(dirL, header=F,colClasses="character")
  
  if (anyDuplicated(mL[,1])!=0)
    stop("Duplicated metabolite list.\n")
  
  checkmetaList <- which(substr(mL[,1],1,1) != "%")
  matchmL <- which(mL[checkmetaList,1] %in% (b[,1]))
  if (length(checkmetaList) != length(matchmL))
  {
    cat("No match found in multi_data(_user).csv for the following metabolite(s):\n")
    print(mL[checkmetaList[setdiff(1:25, matchmL)],1])
    stop("Possible error in file 'metabolitesList.csv'.\n")
  }
    
  
  mL<-mL[,1,drop=FALSE]
  
  write.table(mL,file=dir3, sep = "\t", row.names=FALSE,col.names=FALSE,quote=FALSE)
  ## get spectra data
  if (!missing(BrukerDataDir))
  {
    sa<-readBruker(BrukerDataDir)  
    write.table(sa,file=dir2,row.names=FALSE,col.names=TRUE,quote=FALSE,sep = "\t")
  } else if (!missing(txtFile)) {
    sa<-read.table(txtFile, header=TRUE,sep="\t")
    write.table(sa,file=dir2,row.names=FALSE,col.names=TRUE,quote=FALSE,sep = "\t")
  } else if (!missing(rData)) {
    sa<-rData
    write.table(sa,file=dir2,row.names=FALSE,col.names=TRUE,quote=FALSE,sep = "\t")
  } else {
    sa<-read.table(dir2, header=TRUE,sep="\t")
  }
  
  if ((ncol(sa)-1)<length(sno))
    return(cat("No. of spectra included smaller than input spectra.\n"))
  if (!is.null(colnames(sa))) {
    saname<-colnames(sa)
    saname<-saname[2:length(saname)]
    stit<-rbind(1:length(saname),saname)
  } else {
    stit<-rbind(1:(ncol(sa)-1),1:(ncol(sa)-1))
  }
  ## spectra number
  write.table(stit[,sno],file=paste(dirctime, "/spectraTitle.txt", sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep = "\t")
  
  ## new added, chemshift for each spectrum
  if (chemshif == 1)
  {
    dirCS<-paste(dirA[2],"/chemShiftPerSpec.csv",sep="")
    if (!file.exists(dirCS))
    {
      createChemShiftPerSpec(tempOption = opt, dirA[2])
      stop("No chemShiftPerSpec.csv file found in BatmanInput folder.
				 Creating one now, please modify the values.\n")
		} else {
		  # prepare chemshift for c++
		  chemlist<-read.csv(dirCS,header=T,stringsAsFactors=FALSE,colClasses="character")
		  if (dim(chemlist)[1] != dim(bor)[1])
		    stop(paste("Different number of multiplets in multiplet template file(s) and ", dirCS,
		               "\nPlease modify chemShiftPerSpec.csv or call createChemShiftPerSpec() to create a new one.",sep = ""))
		  for (i in 3:dim(chemlist)[2])
		  {
		    check<-chemlist[,i]=='n'
		    chemlist[check,i]=-50
		  }
		  Dc<-order(chemlist[,1])
		  chemlistor<-chemlist[Dc,]
		  write.table(chemlistor,file=dir6,sep = "\t",row.names = FALSE,col.names = FALSE,quote=FALSE)
		  cp<-file.copy(dir6,dir5)
		}
  }
  
  ## choose whether to parallelize spectra
  if (length(sno)>1 && fixeff == 0)    
  {
    cat("\nHow many parallel processes (multicores) do you want to run the multi-spectra analysis?")
    cat("\n(Enter 1 for running them sequentially.)\n")
    cat("\n Parallel processing of multi spectra currently cannot display the progress\n")
    cat(" bar (or any words), if you input is > 1, please be patient for the results :)\n\n")
    wr<- getinput(lowlim=1,highlim=-1)  
  } else {
    wr<-1
  }
  
  cp<-file.copy(dir2,dir5)
  cp<-file.copy(dir4,dir5)
  
  rr<-0
  if (wr>1) {
    cl<-makeCluster(wr, type = "SOCK")
    registerDoSNOW(cl)  
  } else {
    cat ("Percentage completed...\n")
  }
  ## calling c++ for MCMC
  if (wr>1) {              
    stime <- system.time ({
      multispec<- foreach (nospec = sno,.packages="batman") %dopar% {
        pBar <- txtProgressBar(min =0, max = ito, style = 3)
        out<-.Call("batman",filedir,as.integer(bn),as.integer(ito),as.integer(rr),as.integer(nospec-1), pBar, PACKAGE = "batman")
        close( pBar )
      }})[3]} else if (wr == 1 && fixeff == 1) { 
        stime <- system.time ({
          multispec<- for (nospec in 1:1) {
            pBar <- txtProgressBar(min =0, max = ito, style = 3)
            out<-.Call("batman",filedir,as.integer(bn),as.integer(ito),as.integer(rr),as.integer(length(sno)-1),pBar,PACKAGE = "batman")
            close( pBar )
          }})[3]} else { 
            stime <- system.time ({
              multispec<- for (nospec in sno) {
                pBar <- txtProgressBar(min =0, max = ito, style = 3)
                out<-.Call("batman",filedir,as.integer(bn),as.integer(ito),as.integer(rr),as.integer(nospec-1),pBar,PACKAGE = "batman")
                close( pBar )
              }})[3]}   
  cat ("time ")
  print(stime)
  cat (" second.\n")
  if (wr>1) {stopCluster(cl)}
  
  cat("Reading in saved data in folder\n")
  ## read in batman result files
  BM<-readBatmanOutput(dirctime, dirA[2]) 
  cat(BM$outputDir)
  ## plot results
  if (figBatmanFit)
    plotBatmanFit(BM, saveFigDir = dirctime, listMeta = listMeta)
  if (figRelCon) 
    plotRelCon(BM, saveFigDir = dirctime)
  if (figMetaFit)
    plotMetaFit(BM, saveFigDir = dirctime)
  
  file.remove(dir3)
  file.remove(dir4)
  if (chemshif == 1)
  {
    if (file.exists(dir6))
      file.remove(dir6)
  }
  cat("\nCompleted.\n")
  warnRead<-options(warn = warnDef)
  return(BM)
}