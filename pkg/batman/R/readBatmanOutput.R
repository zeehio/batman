readBatmanOutput<-function(dirOP, dirIP,readMetaIndFitSam = TRUE,
                           readMetaTempHR = TRUE, readMetaTemp = TRUE)
{
  ## written by Dr. Jie Hao, Imperial College London
  warnDef<-options("warn")$warn
  warnRead<-options(warn = -1)
  ## reads in batman output data files 
  dirBO <- paste(dirOP,"/batmanOptions.txt",sep="")
  dirsT <- paste(dirOP,"/spectraTitle.txt",sep="")
  dirOPList <- paste(dirOP,"/metabolitesListUsed.txt",sep="")
  if (!file.exists(dirBO))
    return(cat(paste("\nFile ", dirBO, " is missing in input folder ", dirOP, ".\n.")))
  if (!file.exists(dirsT))
    return(cat(paste("\nFile ", dirsT, " is missing in input folder ", dirOP, ".\n.")))
  if (!file.exists(dirOPList))
    return(cat(paste("\nFile ", dirOPList, " is missing in input folder ", dirOP, ".\n.")))
  
  con  <- file(dirBO, open = "r")
  nlines <-dim(read.table(dirBO,sep="\n",comment.char = ""))[1]
  oneLine <- readLines(con, n = nlines, warn = FALSE)
  fL<-substr(oneLine,1,1)
  nL<-which(is.na(match(fL,"%")))
  myVector <- strsplit(oneLine[nL[2]], ":")
  sno <-getSpectraRange(myVector)
  close(con)
  NoSpectra <- length(sno)
  
  specTitle<-read.table(dirsT, header=FALSE,sep="\t",comment.char ="")
  
  if (file.info(dirOPList)$size == 0 )
    return(cat(paste("\nFile", dirOPList, " is empty. Not fitting anything in ppm range.\n." )))
  mL<-read.csv(dirOPList, header=F,colClasses="character")
  mL<-mL[,1,drop=FALSE]
  
  ## reads in batman() output
  r<-NULL
  rH<-NULL
  bet<-NULL
  L<-NULL
  LH<-NULL
  the<-NULL
  lam<-NULL
  del<-NULL
  delmean<-NULL
  betasam<-NULL
  sfitsam<-NULL
  Ndata <-NULL
  metaindfitsam<-NULL
  jsno <- sno[1]
  fdir <- paste(dirOP, "/specFit_",jsno,"_rr_0.txt", sep="")
  if (file.exists(fdir))
  {                
    if (!file.info(fdir)$size == 0) 
    {
      r<-read.table(fdir,sep = "\t",header=T,comment.char = "", colClasses="numeric")
    } 
  }
  fdir <- paste(dirOP, "/NMRdata_mod_",jsno,".txt", sep="")
  if (file.exists(fdir))
  {                
    if (!file.info(fdir)$size == 0) 
    {
      Ndata<-read.table(fdir, header = FALSE, sep = "\t",comment.char = "", colClasses="numeric")
    }   
  } 
  fdir <- paste(dirOP, "/specFitHR_",jsno,"_rr_0.txt", sep="")
  if (file.exists(fdir))
  {                
    if (!file.info(fdir)$size == 0) 
    {
      rH<-read.table(fdir,sep = "\t",header=F,comment.char = "", colClasses="numeric")
      rH<-cbind(Ndata,rH)
    }
  }
  if (readMetaTempHR)
  {
  fdir <- paste(dirOP, "/metaTempHR_",jsno,"_rr_0.txt", sep="")
  if (file.exists(fdir))
  {                
    if (!file.info(fdir)$size == 0) 
    {
      LH<-read.table(fdir,sep = "\t", header=F,comment.char = "", colClasses="numeric")
      if (length(mL) == ncol(LH))
        names(LH)<-mL
      if (nrow(mL) == ncol(LH))
        names(LH)<-t(mL)
    }
  }
  }
  if (readMetaTemp)
  {
  fdir <- paste(dirOP, "/metaTemp_",jsno,"_rr_0.txt", sep="")
  if (file.exists(fdir))
  {                
    if (!file.info(fdir)$size == 0) 
    {
      L<-read.table(fdir,sep = "\t", header=F,comment.char = "")
      if (length(mL) == ncol(L))
        names(L)<-mL
      if (nrow(mL) == ncol(L))
        names(L)<-t(mL)
    }
  }
  }
  fdir <- paste(dirOP, "/beta_",jsno,"_rr_0.txt", sep="")
  if (file.exists(fdir))
  {                
    if (!file.info(fdir)$size == 0) 
    {
      bet<-read.table(fdir,header=F,comment.char = "")
      if (length(mL) == nrow(bet))
        row.names(bet)<-mL
      if ( nrow(mL) == nrow(bet) )
        row.names(bet)<-t(mL)
    }
  }
  fdir <- paste(dirOP, "/theta_sam_",jsno,"_rr_0.txt", sep="") 
  if (file.exists(fdir))
  {                
    if (!file.info(fdir)$size == 0) 
    {
      if (!(file.info(fdir)$size == 0))
      {
        the<-read.table(fdir,sep = "\t",header=F,comment.char = "") 
      }
    }
  }
  fdir <- paste(dirOP, "/metaFit_sam_",jsno,"_rr_0.txt", sep="") 
  if (file.exists(fdir))
  {                
    if (!file.info(fdir)$size == 0) 
    {
      sfitsam<-read.table(fdir,sep = "\t",header=F,comment.char = "")
    }  
  } 
  fdir <- paste(dirOP, "/lambda_sam_",jsno,"_rr_0.txt", sep="")
  if (file.exists(fdir))
  {                
    if (!file.info(fdir)$size == 0) 
    {
      lam<-read.table(fdir,header=F,comment.char = "")
    }    
  }   
  fdir <- paste(dirOP, "/beta_sam_",jsno,"_rr_0.txt", sep="")
  if (file.exists(fdir))
  {                
    if (!file.info(fdir)$size == 0) 
    {
      betasam<-read.table(fdir,sep = "\t",header=F,comment.char = "", colClasses="numeric")

      if (length(mL) == nrow(betasam))
        row.names(betasam)<-mL
      if (nrow(mL) == nrow(betasam))
        row.names(betasam)<-t(mL)
    }  
  }
  fdir <- paste(dirOP, "/delta_sam_",jsno,"_rr_0.txt", sep="")
  if (file.exists(fdir))
  {             
    if (!file.info(fdir)$size == 0) 
    {
      del<-read.delim(fdir,header=T)
    }
  }
  fdir <- paste(dirOP, "/delta_sam_",jsno,".txt", sep="")
  if (file.exists(fdir))
  {             
    if (!file.info(fdir)$size == 0) 
    {
      del<-read.delim(fdir,header=T)
    }
  }
  fdir <- paste(dirOP, "/delta_draw_mean_",jsno,".txt", sep="")
  if (file.exists(fdir))
  {                
    if (!file.info(fdir)$size == 0) 
    {
      delmean<-read.delim(fdir,header=F)
      if (ncol(del) == nrow(delmean))
        row.names(delmean)<-names(del)
    }
  }
  if (readMetaIndFitSam)
  {
  ## read in individual metabolite fit posteriors
  brow<-nrow(betasam)
  bcol<-ncol(betasam)
  fdir <- paste(dirOP, "/metaIndFit_sam_",jsno,"_rr_0.txt", sep="")
  if (file.exists(fdir))
  { 
    if (!file.info(fdir)$size == 0)
    {
      metafitsamInd<-read.table(fdir,sep="\t",header=F,comment.char = "", colClasses="numeric")  
      for (j in 1:bcol)
      {
        if (j == 1)
          metaindfitsam<-t(t(metafitsamInd[,((j-1)*brow+1):(j*brow)])*betasam[,j])
        else
          metaindfitsam<-cbind(metaindfitsam, t(t(metafitsamInd[,((j-1)*brow+1):(j*brow)])*betasam[,j]))
      }                  
    }
  }
  }
  ## more than 1 spectra
  if (NoSpectra>1)
  {
    for (i in sno[2:NoSpectra])   
    { 
      fdir <- paste(dirOP, "/specFit_", i,"_rr_0.txt", sep="")
      if (file.exists(fdir))
      {
        if (!file.info(fdir)$size == 0)
        {
          r<-cbind(r, read.table(fdir,sep = "\t",header=T,comment.char = ""))
        }
      }    
      fdir <- paste(dirOP, "/NMRdata_mod_",i,".txt", sep="")
      if (file.exists(fdir))
      {
        if (!file.info(fdir)$size == 0)
        {
          Ndata<-read.table(fdir, header = FALSE, sep = "\t",comment.char = "")
        }
      }
      fdir <- paste(dirOP, "/specFitHR_", i,"_rr_0.txt", sep="")
      if (file.exists(fdir))
      { 
        if (!file.info(fdir)$size == 0)
        {               
          rH<-cbind(rH, Ndata, read.table(fdir,sep = "\t",header=F,comment.char = ""))
        }
      }
      if (readMetaTempHR)
      {
      fdir <- paste(dirOP,  "/metaTempHR_",i, "_rr_0.txt", sep="")
      if (file.exists(fdir))
      {
        if (!file.info(fdir)$size == 0)
        {
          LtH<-read.table(fdir,sep = "\t",header=F,comment.char = "", colClasses="numeric")
          if (length(mL) == ncol(LtH))
            names(LtH)<-mL 
          if (nrow(mL) == ncol(LtH))
            names(LtH)<-t(mL) 
          LH<-cbind(LH, LtH)
        }
      }
      }
      fdir <- paste(dirOP,  "/beta_",i, "_rr_0.txt", sep="")
      if (file.exists(fdir))
      {
        if (!file.info(fdir)$size == 0)
        {
          bet<-cbind(bet, read.table(fdir,header=F,comment.char = "", colClasses="numeric"))
        }
      }
      if (readMetaTemp)
      {
      fdir <- paste(dirOP,  "/metaTemp_",i, "_rr_0.txt", sep="")
      if (file.exists(fdir))
      {
        if (!file.info(fdir)$size == 0)
        {
          Lt<-read.table(fdir,sep = "\t",header=F,comment.char = "", colClasses="numeric")
          if (length(mL) == ncol(Lt))
            names(Lt)<-mL 
          if (nrow(mL) == ncol(Lt))
            names(Lt)<-t(mL) 
          L<-cbind(L, Lt)
        }
      }
      }
      fdir <- paste(dirOP, "/beta_sam_",i, "_rr_0.txt", sep="")
      if (file.exists(fdir))
      {
        if (!file.info(fdir)$size == 0) 
        {
          betasamInd<-read.table(fdir,sep = "\t", header=F,comment.char = "", colClasses="numeric")
          betasam<-cbind(betasam, betasamInd)
        }
      }                
      fdir <- paste(dirOP, "/delta_draw_mean_",i,".txt", sep="")
      if (file.exists(fdir))
      {
        if (!file.info(fdir)$size == 0) 
        {
          delmean<-cbind(delmean, read.delim(fdir,header=F))
        }
      }
      fdir <- paste(dirOP,  "/delta_sam_",i, "_rr_0.txt", sep="")
      if (file.exists(fdir))
      {
        if (!file.info(fdir)$size == 0) 
        {
          del<-cbind(del, read.delim(fdir, header=T))
        }
      }
      fdir <- paste(dirOP,  "/delta_sam_",i, ".txt", sep="")
      if (file.exists(fdir))
      {
        if (!file.info(fdir)$size == 0) 
        {
          del<-cbind(del, read.delim(fdir, header=T))
        }
      }
      fdir <- paste(dirOP,  "/theta_sam_",i, "_rr_0.txt", sep="")
      if (file.exists(fdir))
      {                if (!file.info(fdir)$size == 0) 
      {
        the<-cbind(the, read.table(fdir, sep = "\t", header=F,comment.char = "", colClasses="numeric"))
      } 
      }
      fdir <- paste(dirOP,  "/metaFit_sam_",i, "_rr_0.txt", sep="")
      if (file.exists(fdir))
      {    
        if (!file.info(fdir)$size == 0) 
        {
          sfitsam<-cbind(sfitsam, read.table(fdir, sep = "\t", header=F,comment.char = "", colClasses="numeric"))
        }
      }
      fdir <- paste(dirOP, "/lambda_sam_",i,"_rr_0.txt", sep="")
      if (file.exists(fdir))
      {    
        if (!file.info(fdir)$size == 0) 
        {
          lam<-cbind(lam,read.table(fdir,header=F,comment.char = "", colClasses="numeric"))
        }  
      }
      if (readMetaIndFitSam)
      {
      fdir <- paste(dirOP, "/metaIndFit_sam_",i,"_rr_0.txt", sep="") 
      if (file.exists(fdir))
      {
        if (!file.info(fdir)$size == 0)
        {
          metafitsamInd<-read.table(fdir,sep="\t",header=F,comment.char = "", colClasses="numeric") 
          for (j in 1:bcol)
            metaindfitsam<-cbind(metaindfitsam, t(t(metafitsamInd[,((j-1)*brow+1):(j*brow)])*betasamInd[,j]))                
        } 
      }
      }
    }
  }
  ## column name for beta
  if (!is.null(bet))
  {
    if (ncol(bet) == ncol(specTitle))
      names(bet)<- t(specTitle[2,])
    else
      names(bet)<- t(specTitle[2,1:dim(bet)[2]])
  }
  ## column name for delta
  if (!is.null(delmean))
  {
    if (ncol(delmean) == ncol(specTitle))
      names(delmean)<- t(specTitle[2,])
    else
      names(delmean)<- t(specTitle[2,1:dim(delmean)[2]])
  }
  ## reads in batmanrerun() output
  rrr<-NULL
  rrrH<-NULL
  betrr<-NULL
  Lrr<-NULL
  LrrH<-NULL
  betasamrr<-NULL
  therr<-NULL
  sfitsamrr<-NULL
  metaindfitsamrr<-NULL
  rr <- 1
  
  fdir <- paste(dirOP, "/specFit_",jsno,"_rr_1.txt", sep="")
  if (file.exists(fdir))
  {                
    if (!file.info(fdir)$size == 0) 
    {
      rrr<-read.table(fdir,sep = "\t",header=T,comment.char = "", colClasses="numeric")
    }
  }
  fdir <- paste(dirOP, "/NMRdata_mod_1.txt", sep="")
  if (file.exists(fdir))
  {                
    if (!file.info(fdir)$size == 0) 
    {
      Ndata<-read.table(fdir, header = FALSE, sep = "\t",comment.char = "", colClasses="numeric")
    }
  }
  fdir <- paste(dirOP, "/specFitHR_",jsno,"_rr_",rr,".txt", sep="")
  if (file.exists(fdir))
  {                
    if (!file.info(fdir)$size == 0) 
    {
      rrrH<-read.table(fdir,sep = "\t",header=F,comment.char = "", colClasses="numeric")
      rrrH<-cbind(Ndata,rrrH)
    }
  }
  if (readMetaTempHR)
  {
  fdir <- paste(dirOP, "/metaTempHR_",jsno,"_rr_",rr,".txt", sep="")
  if (file.exists(fdir))
  {                
    if (!file.info(fdir)$size == 0) 
    {
      LrrH<-read.table(fdir,sep = "\t",header=F,comment.char = "", colClasses="numeric")
      if (length(mL) == ncol(LrrH))
        names(LrrH)<-mL
      if ( nrow(mL) == ncol(LrrH))
        names(LrrH)<-t(mL)
    }
  }
  }
  if (readMetaTemp)
  {
  fdir <- paste(dirOP, "/metaTemp_",jsno,"_rr_",rr,".txt", sep="")
  if (file.exists(fdir))
  {                
    if (!file.info(fdir)$size == 0) 
    {
      Lrr<-read.table(fdir,sep = "\t",header=F,comment.char = "", colClasses="numeric")
      if (length(mL) == ncol(Lrr))
        names(Lrr)<-mL
      if ( nrow(mL) == ncol(Lrr))
        names(Lrr)<-t(mL)
    }
  }
  }
  fdir <- paste(dirOP, "/beta_",jsno,"_rr_",rr,".txt", sep="")
  if (file.exists(fdir)) 
  {                
    if (!file.info(fdir)$size == 0) 
    {
      betrr<-read.table(fdir,header=F,comment.char = "", colClasses="numeric")
      if (length(mL) == nrow(betrr))
        row.names(betrr)<-mL
      if (nrow(mL) == ncol(betrr))
        row.names(betrr)<-t(mL)
    }
  }
  fdir <- paste(dirOP, "/beta_sam_",jsno,"_rr_",rr,".txt", sep="")
  if (file.exists(fdir)) 
  {
    if (!file.info(fdir)$size == 0) 
    {
      betasamrr<-read.table(fdir,sep = "\t",header=F,comment.char = "", colClasses="numeric")
      if (length(mL) == nrow(betasamrr))
        row.names(betasamrr)<-mL
      if (nrow(mL) == nrow(betasamrr))
        row.names(betasamrr)<-t(mL)
    }
  }
  fdir <- paste(dirOP, "/theta_sam_",jsno,"_rr_",rr,".txt", sep="")
  if (file.exists(fdir))
  {
    if (!file.info(fdir)$size == 0) 
    {
      therr<-read.table(fdir,sep = "\t",header=F,comment.char = "", colClasses="numeric")
    }
  }
  fdir <- paste(dirOP, "/metaFit_sam_",jsno,"_rr_",rr,".txt", sep="")
  if (file.exists(fdir))
  {if (!file.info(fdir)$size == 0) 
  {
    sfitsamrr<-read.table(fdir,sep = "\t",header=F,comment.char = "", colClasses="numeric")
  }}
  if (readMetaIndFitSam)
  {
  ## individual metabolite fit posteriors
  brow<-nrow(betasamrr)
  bcol<-ncol(betasamrr)
  fdir <- paste(dirOP, "/metaIndFit_sam_",jsno,"_rr_",rr,".txt", sep="")  
  if (file.exists(fdir))
  { 
    if (!file.info(fdir)$size == 0)
    {
      metafitsamrrInd<-read.table(fdir,sep="\t",header=F,comment.char = "", colClasses="numeric") 
      for (j in 1:bcol)
      {
        if (j == 1)
          metaindfitsamrr<-t(t(metafitsamrrInd[,((j-1)*brow+1):(j*brow)])*betasamrr[,j])
        else
          metaindfitsamrr<-cbind(metaindfitsamrr, t(t(metafitsamrrInd[,((j-1)*brow+1):(j*brow)])*betasamrr[,j]))
      }                       
    }
  }
  }
  ## more spectra
  if (NoSpectra>1)
  {
    for (i in sno[2:NoSpectra])   
    { 
      fdir <- paste(dirOP, "/specFit_", i,"_rr_",rr,".txt", sep="")
      if (file.exists(fdir))
      {
        if (!file.info(fdir)$size == 0) 
        {
          rrr<-cbind(rrr, read.table(fdir,sep = "\t",header=T,comment.char = "", colClasses="numeric"))
        }
      }           
      fdir <- paste(dirOP, "/NMRdata_mod_",i,".txt", sep="")
      if (file.exists(fdir))
      {
        if (!file.info(fdir)$size == 0) 
        {
          Ndata<-read.table(fdir, header = FALSE, sep = "\t",comment.char = "", colClasses="numeric")
        }
      }
      fdir <- paste(dirOP, "/specFitHR_", i,"_rr_",rr,".txt", sep="")
      if (file.exists(fdir))
      {
        if (!file.info(fdir)$size == 0) 
        {
          rrrH<-cbind(rrrH, Ndata, read.table(fdir,sep = "\t",header=F,comment.char = "", colClasses="numeric"))
        }
      }
      if (readMetaTempHR)
      {
      fdir <- paste(dirOP, "/metaTempHR_",i, "_rr_",rr,".txt", sep="")
      if (file.exists(fdir))
      {
        if (!file.info(fdir)$size == 0) 
        {
          LrrtH<-read.table(fdir,sep = "\t",header=F,comment.char = "", colClasses="numeric")
          if (length(mL) == ncol(LrrtH))
            names(LrrtH)<-mL
          if ( nrow(mL) == ncol(LrrtH))
            names(LrrtH)<-t(mL)
          LrrH<-cbind(LrrH, LrrtH)
        }
      }
      }
      if (readMetaTemp)
      {
      fdir <- paste(dirOP, "/metaTemp_",i, "_rr_",rr,".txt", sep="")
      if (file.exists(fdir))
      {
        if (!file.info(fdir)$size == 0) 
        {
          Lrrt<-read.table(fdir,sep = "\t",header=F,comment.char = "", colClasses="numeric")
          if (length(mL) == ncol(Lrrt))
            names(Lrrt)<-mL
          if ( nrow(mL) == ncol(Lrrt))
            names(Lrrt)<-t(mL)
          Lrr<-cbind(Lrr, Lrrt)
        }   
      }
      }
      fdir <- paste(dirOP, "/beta_",i, "_rr_",rr,".txt", sep="")
      if (file.exists(fdir)) 
      {
        if (!file.info(fdir)$size == 0) 
        {
          betrr<-cbind(betrr, read.table(fdir,header=F,comment.char = "", colClasses="numeric"))
        }
      }            
      fdir <- paste(dirOP, "/beta_sam_",i, "_rr_",rr,".txt", sep="")
      if (file.exists(fdir)) 
      {
        if (!file.info(fdir)$size == 0) 
        {
          betasamrrInd<-read.table(fdir,sep = "\t",header=F,comment.char = "", colClasses="numeric")
          betasamrr<-cbind(betasamrr, betasamrrInd)  
        }
      }
      fdir <- paste(dirOP,  "/theta_sam_",i, "_rr_",rr,".txt", sep="")
      if (file.exists(fdir))
      {
        if (!file.info(fdir)$size == 0) 
        {
          therr<-cbind(therr, read.table(fdir, sep = "\t", header=F,comment.char = "", colClasses="numeric"))    
        }
      }
      fdir <- paste(dirOP,  "/metaFit_sam_",i, "_rr_",rr,".txt", sep="")
      if (file.exists(fdir))
      {
        if (!file.info(fdir)$size == 0) 
        {
          sfitsamrr<-cbind(sfitsamrr, read.table(fdir, sep = "\t", header=F,comment.char = "", colClasses="numeric"))    
        }   
      }       
      if (readMetaIndFitSam)
      {
      fdir <- paste(dirOP, "/metaIndFit_sam_",i, "_rr_",rr,".txt", sep="") 
      if (file.exists(fdir))
      {
        if (!file.info(fdir)$size == 0) 
        {
          metafitsamrrInd<-read.table(fdir,sep="\t",header=F,comment.char = "", colClasses="numeric")   
          for (j in 1:bcol)
            metaindfitsamrr<-cbind(metaindfitsamrr, t(t(metafitsamrrInd[,((j-1)*brow+1):(j*brow)])*betasamrrInd[,j]))                    
        }
      }
      }
    }    
  }
  ## column name for beta
  if (!is.null(betrr))
  {
    if (ncol(betrr) == ncol(specTitle))
      names(betrr)<- t(specTitle[2,])
    else 
      names(betrr)<- t(specTitle[2,1:dim(betrr)[2]])
    row.names(betrr)<-row.names(bet)
    
  }
  if (!is.null(bet))
  {
    bet2<-bet
    Metabolite<-row.names(bet2)
    row.names(bet2)<-NULL
    bet3<-cbind(Metabolite,bet2)
    dirOPRC <- paste(dirOP,"/RelCon.txt",sep="")
    if (!file.exists(dirOPRC))
      write.table(bet3,file = dirOPRC, sep = "\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
  }
  if (!is.null(betrr))
  {
    betrr2<-betrr
    Metabolite<-row.names(betrr2)
    row.names(betrr2)<-NULL
    betrr3<-cbind(Metabolite,betrr2)
    dirOPRCr<- paste(dirOP,"/RelConRerun.txt",sep="")
    if (!file.exists(dirOPRCr))
      write.table(betrr3,file = dirOPRCr,sep = "\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
  }
  if (!is.null(delmean))
  {
    delmean2<-delmean
    Multiplet<-row.names(delmean2)
    row.names(delmean2)<-NULL
    delmean3<-cbind(Multiplet,delmean2)
    dirOPMP<- paste(dirOP,"/MultipletsPpmShifts.txt",sep="")
    if (!file.exists(dirOPMP))
      write.table(delmean3,file = dirOPMP, sep = "\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
  }
  ns<-5
  betvA <-NULL  
  betvArr <-NULL 
  vA<-NULL
  vArr<-NULL
  dirOPRCCI<- paste(dirOP,"/RelConCreInt.txt",sep="")
  if (!file.exists(dirOPRCCI))
  {
    if (!is.null(betasam))
    {
      f<-nrow(betasam)
      fc<-ncol(betasam)
      #sno<-length(r)/ns
      #ind<-fc/length(sno)
      ind<-fc/dim(bet)[2]
      for (i in 1:dim(bet)[2]) 
      {
        vAll<-NULL
        for (j in 1:f)
        {
          v<-quantile(betasam[j,((i-1)*ind+1):(i*ind)],p = c(2.5,97.5)/100)  
          vAll<-rbind(vAll,v)
        }
        if (i == 1)
          vA<-vAll
        else
          vA<-cbind(vA,vAll)
      }   
      vA2<-vA
      percentage<-names(vA2)
      Metabolite<-rep(t(specTitle[2,]),each=2)
      vA2<-rbind(Metabolite,percentage,vA2)
      row.names(vA2)<-c("Metabolite", "Percentage",row.names(bet))
      write.table(vA2,file = dirOPRCCI, sep = "\t",row.names = TRUE,col.names = FALSE,quote=FALSE)
    }
  }
  else
  {
    vA<- read.table(file = dirOPRCCI, sep = "\t",header = TRUE,comment.char = "")
    vA2 <- vA[2:dim(vA)[1], 2:dim(vA)[2]]
    row.names(vA2)<-vA[2:dim(vA)[1],1]
    names(vA2)<-t(vA[1,2:dim(vA)[2]])
    vA <- vA2
  }
  dirOPRCCIr<- paste(dirOP,"/RelConCreIntRerun.txt",sep="")
  if (!file.exists(dirOPRCCIr))
  {
    if (!is.null(betasamrr))
    { 
      f<-nrow(betasamrr)
      fc<-ncol(betasamrr)
      #ind<-fc/length(sno)
      ind<-fc/dim(betrr)[2]
      for (i in 1:dim(betrr)[2]) 
      {
        vAllrr<-NULL
        for (j in 1:f)
        {
          v<-quantile(betasamrr[j,((i-1)*ind+1):(i*ind)],p = c(2.5,97.5)/100)  
          vAllrr<-rbind(vAllrr,v)
        }
        if (i == 1)
          vArr<-vAllrr
        else
          vArr<-cbind(vArr,vAllrr)
      }   
      vArr2<-vArr
      percentage<-names(vArr2)
      Metabolite<-rep(t(specTitle[2,]),each=2)
      vArr2<-rbind(Metabolite,percentage,vArr2)
      row.names(vArr2)<-c("Metabolite", "Percentage",row.names(bet))
      write.table(vArr2, file = dirOPRCCIr, sep = "\t", row.names = TRUE, col.names = FALSE, quote=FALSE)
    }  
  }
  else
  {
    vArr<- read.table(file = dirOPRCCIr, sep = "\t",header = TRUE,comment.char = "")
    vArr2 <- vArr[2:dim(vArr)[1], 2:dim(vArr)[2]]
    row.names(vArr2)<-vArr[2:dim(vArr)[1],1]
    names(vArr2)<-t(vArr[1,2:dim(vArr)[2]])
    vArr <- vArr2
  }
  if (missing(dirIP))
    dirIP <- NULL
  ## batman fitting results
  BM<-list(specTitle = specTitle, specRange = sno, sFit=r,sFitHR=rH, beta=bet, betaSam=betasam, betaCI=vA, metaTemp=L,
           metaTempHR=LH, metaFitSam = sfitsam, metaIndFitSam = metaindfitsam, thetaSam=the, delta = delmean, deltaSam = del,
           sFitRerun=rrr, sFitRerunHR=rrrH, betaRerun = betrr, betaSamRerun = betasamrr, betaCIRerun=vArr,
           metaTempRerun = Lrr, metaTempRerunHR = LrrH, metaFitSamRerun=sfitsamrr, metaIndFitSamRerun=metaindfitsamrr,
           thetaSamRerun=therr,outputDir = dirOP, inputDir = dirIP) 
  warnRead<-options(warn = warnDef)
  return (BM)
}    

