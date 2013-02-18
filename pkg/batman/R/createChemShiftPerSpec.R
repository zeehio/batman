createChemShiftPerSpec <-function(templateOption = 1, dirIP)
{
## written by Dr. Jie Hao, Imperial College London
  if(missing(dirIP))
  {
    return("Please provide input to 'dirIP'.\n")
  }
  dirA <- dirIP
  
  dirCS<-paste(dirA,"/chemShiftPerSpec.csv",sep="")
  dir1<-paste(dirA,"/batmanOptions.txt",sep="")
  dir2<-paste(dirA,"/NMRdata.txt",sep="")
  #dirL<-paste(dirA,"/metabolitesList.csv",sep="")
  dirR<-paste(dirA,"/multi_data.csv",sep="")
  dirRU<-paste(dirA,"/multi_data_user.csv",sep="")
  
    if (templateOption == 1)
    {
      cat("Copying multiplet list from multi_data.csv to chemShiftPerSpec.csv...\n")
      b<-read.csv(dirR,header=T,stringsAsFactors=FALSE,colClasses="character")
    } else if (templateOption == 2)  {
      cat("Copying multiplet list from multi_data_user.csv to chemShiftPerSpec.csv...\n")
      b<-read.csv(dirRU,header=T,stringsAsFactors=FALSE,colClasses="character")
    } else if (templateOption == 3)  {
      cat("Copying multiplet list from multi_data.csv and multi_data_user.csv to chemShiftPerSpec.csv...\n")
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
     ##D<-order(b[,1])
    tempData<-b
  
  ## read in batman optitons
  con  <- file(dir1, open = "r")
  oneLine <- readLines(con, n = 30, warn = FALSE)
  fL<-substr(oneLine,1,1)
  nL<-which(is.na(match(fL,"%")))
  myVector <- strsplit(oneLine[nL[2]], ":")
	ranges <- strsplit(myVector[[1]][2],",")
	sno <- NULL
	sno <-getSpectraRange(myVector)
##  sNo <- as.numeric(myVector[[1]][2])
  close(con)
  sNo <- length(sno)

	sa<-read.table(dir2, header=TRUE,sep="\t",comment.char = "")
    if ((ncol(sa)-1)<sNo)
      return(cat("No. of spectra included smaller than input spectra.\n"))
    if (!is.null(colnames(sa))) {
      saname<-colnames(sa)
      saname<-saname[2:length(saname)]
      sTitle<-rbind(1:length(saname),saname)
    } else {
      sTitle<-rbind(1:(ncol(sa)-1),1:(ncol(sa)-1))
    }
  
  
  chemlist <- NULL
  #for (n in 1:length(metaList))
  #{
    #chemlist<-rbind(chemlist,tempData[which(tempData[,1] == metaList[n,]),1:2])
  #}
  chemlist<-tempData[,1:2]
  #dim(chemlist)<-c(length(chemlist),1)
  nlist <- matrix("n", dim(chemlist)[1], 1)
  for (n in 1:dim(sTitle)[2])
  {
    chemlist<-cbind(chemlist, nlist)
  }
  #return(list(chemlist = chemlist, sTitle = sTitle[2,]))
  names(chemlist)<-c("multiplets","pos_in_ppm", sTitle[2,])
  write.table(chemlist,file=dirCS,sep = ",",row.names = FALSE,col.names = TRUE)
  
}