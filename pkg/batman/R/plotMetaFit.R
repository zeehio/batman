plotMetaFit<-function(BM, from, to, metaName, saveFig = TRUE,saveFigDir = BM$outputDir, 
                      prefixFig, rerun = FALSE, overwriteFig = FALSE)
{
  ## written by Dr. Jie Hao
  ## plot metabolites fit posteriors with 95% credible Interval
  if (missing(BM))
    return(cat("Please input batman data list.\n"))
  
  warnDef<-options("warn")$warn
  warnRead<-options(warn = -1)
  ptype = "pdf"
  n <- 2
  ns<-5
  
  if (missing(from))
    from = min(BM$sFit[1,1],BM$sFit[nrow(BM$sFit),1])
  if (missing(to))
    to =  max(BM$sFit[1,1],BM$sFit[nrow(BM$sFit),1])
  
  if (from > to)
  {
    temp<-from
    from<-to
    to<-temp
  }
  
  pind<-which(BM$sFit[,1]<=to & BM$sFit[,1]>=from)
  if (length(pind)== 0)
    pind<-which(BM$sFit[,1]<=from & BM$sFit[,1]>=to)
  
  m<-row.names(BM$beta)
  if (!missing(metaName))
  {
    mind<-which(!is.na(match(tolower(m),tolower(metaName)))) 
    if (length(mind)==0) 
      cat("No matching metabolite found...\n")    
  } else {
    metaName<-NULL
    mind<-NULL
  }
  ## plot batman result
  if (!is.null(BM$metaFitSam) && !rerun) 
  {
    f<-nrow(BM$metaFitSam)
    fc<-ncol(BM$metaFitSam)
    
    sno<- BM$specRange
    ind<-fc/length(sno)
    outpdf1<-NULL
    for (i in 1:length(sno)) 
    {
      ## set subplot
      if ((i%%n) == 1){ 
        if ((length(sno)-i)>=(n-1)) {
          if (!missing(prefixFig))
            outpdf1 <- paste(saveFigDir, "/", prefixFig,"_spec_",sno[i],"to",sno[i+n-1],"_mFitSam_", metaName, ".",ptype, sep="")
          else
            outpdf1 <- paste(saveFigDir,"/spec_",sno[i],"to",sno[i+n-1],"_mFitSam_", metaName, ".",ptype, sep="")
          x11(15,7)
          par(mfrow=c(n,1))  
        } else {
          if (!missing(prefixFig))
            outpdf1 <- paste(saveFigDir, "/", prefixFig, "_spec_",sno[i],"_mFitSam_", metaName, ".",ptype, sep="")
          else
            outpdf1 <- paste(saveFigDir,"/spec_",sno[i],"_mFitSam_", metaName, ".",ptype, sep="")
          x11(15,7)
        }
      } 
      v<-NULL
      j = ns*(i-1)+1
      ## set plot parameters
      if (!missing(metaName) && length(mind)!= 0 )
      {
        brow<-nrow(BM$betaSam)
        temp<-NULL
        temp<-BM$metaIndFitSam[,seq(mind,ncol(BM$metaIndFitSam),brow)]
        metaFittemp<-temp[,((i-1)*ind+1):(i*ind)]       
        maintitle <- paste("NMR Spectrum ", sno[i], ": ",BM$specTitle[2,i], "\n", metaName, " fit sample", sep="")  
        metaFitFinal <-(BM$beta[mind,i]*BM$metaTemp[pind,mind+(i-1)*nrow(BM$beta)])
      } else {
        metaFittemp<-BM$metaFitSam[,((i-1)*ind+1):(i*ind)]
        maintitle <- paste("NMR Spectrum ", sno[i], ": ",BM$specTitle[2,i], "\nMetabolite fit sample", sep="")
        metaFitFinal <-BM$sFit[pind,j+2]
      }
      if (length(((i-1)*ind+1):(i*ind))>1)
      {
        v<-apply(metaFittemp, 1, quantile, probs = c(2.5,97.5)/100,  na.rm = TRUE) 
        v<-t(v)
      } else {
        v<-cbind(metaFittemp,metaFittemp)
      }
      ## plot
      plot(BM$sFit[pind,j], BM$sFit[pind,j+1],col="blue", lwd = 3, lty = 1, type = "l", 
           xlim=rev(range(BM$sFit[pind,j])),ylim = c(min(v[pind,1]),max(max(v[pind,2]),BM$sFit[pind,j+1])), 
           xlab ="ppm", ylab = "Standardized Intensity", 
           main = maintitle) 
      #axis(1, 1:length(pind), lab = format(BM$sFit[pind,j]))
      lines(BM$sFit[pind,j], v[pind,1],col = "black", lwd = 3, lty = 4)
      lines(BM$sFit[pind,j], v[pind,2],col = "gray", lwd = 3, lty = 4)
      lines(BM$sFit[pind,j], metaFitFinal,col="green", lwd = 3, lty = 1)
      
      legend("topright", legend = c("Metabolite fit","2.5% quantiles", "97.5% quantiles", "Original Spectrum"),
             col = c("green", "black", "gray", "blue"), pt.cex=3, lty = c(1,4,4,1),lwd = c(3,3,3,3), cex = 0.8)
      ## save
      if ((length(sno) == i || !(i%%n)) && saveFig) {
        if (file.exists(outpdf1) && !overwriteFig)
          cat("Can't save figure, file", outpdf1, "already exists.\n")
        else
          df = dev.copy2pdf(device=x11, file = outpdf1)
      }
    }
  }
  ## plot rerun batman result
  else if (!is.null(BM$metaFitSamRerun) && rerun) 
  {     
    f<-nrow(BM$metaFitSamRerun)
    fc<-ncol(BM$metaFitSamRerun)
    
    sno<-BM$specRange
    ind<-fc/length(sno)
    outpdf2<-NULL
    for (i in 1:length(sno)) {
      ## set subplot
      if ((i%%n) == 1){ 
        if ((length(sno)-i)>=(n-1)) {
          if (!missing(prefixFig))
            outpdf2 <- paste(saveFigDir, "/", prefixFig,"_specRerun_",sno[i],"to",sno[i+n-1],"_mFitSam_", metaName, ".",ptype, sep="")
          else
            outpdf2 <- paste(saveFigDir,"/specRerun_",sno[i],"to",sno[i+n-1],"_mFitSam_", metaName, ".",ptype, sep="")
          x11(15,7)
          par(mfrow=c(n,1))	
        } else {
          if (!missing(prefixFig))
            outpdf2 <- paste(saveFigDir, "/", prefixFig, "_specRerun_",sno[i],"_mFitSam_", metaName,".",ptype, sep="")
          else
            outpdf2 <- paste(saveFigDir,"/specRerun_",sno[i],"_mFitSam_", metaName,".",ptype, sep="")
          x11(15,7)
        }
      }     	
      v<-NULL
      j = ns*(i-1)+1
      ## set plot parameters
      if (!missing(metaName) && length(mind)!= 0 )
      {
        brow<-nrow(BM$betaSamRerun)
        temp<-NULL
        temp<-BM$metaIndFitSamRerun[,seq(mind,ncol(BM$metaIndFitSamRerun),brow)]
        metaFittemp<-temp[,((i-1)*ind+1):(i*ind)]      
        maintitle <- paste("NMR Spectrum ", sno[i], ": ", BM$specTitle[2,i], "(Rerun)\n", metaName, " fit sample", sep="")  
        metaFitFinal <-(BM$betaRerun[mind,i]*BM$metaTempRerun[pind,mind+(i-1)*nrow(BM$betaRerun)])
      } else {
        metaFittemp<-BM$metaFitSamRerun[,((i-1)*ind+1):(i*ind)]
        maintitle <- paste("NMR Spectrum ", sno[i], ": ", BM$specTitle[2,i], "(Rerun)\nMetabolite fit sample", sep="")
        metaFitFinal<-BM$sFitRerun[pind,j+2]
      }
      if (length(((i-1)*ind+1):(i*ind))>1)
      {
        v<-apply(metaFittemp, 1, quantile, probs = c(2.5,97.5)/100,  na.rm = TRUE) 
        v<-t(v)
      } else {
        v<-cbind(metaFittemp,metaFittemp)
      }
      plot(BM$sFitRerun[pind,j], BM$sFitRerun[pind,j+1],col="Blue", lwd = 3, lty = 1, type = "l", 
           xlim=rev(range(BM$sFitRerun[pind,j])),ylim = c(min(v[pind,1]),max(max(v[pind,2]),BM$sFitRerun[pind,j+1])), 
           xlab ="ppm", ylab = "Standardized Intensity", main = maintitle) 
      #axis(1, 1:length(BM$sFitRerun[pind,j]), lab = format(BM$sFitRerun[pind,j]))
      lines(BM$sFitRerun[pind,j], v[pind,1],col = "black", lwd = 3, lty = 4)
      lines(BM$sFitRerun[pind,j], v[pind,2],col = "gray", lwd = 3, lty = 4)
      lines(BM$sFitRerun[pind,j], metaFitFinal,col="green", lwd = 3, lty = 1)
      
      legend("topright", legend = c("Metabolite fit","2.5% quantiles", "97.5% quantiles", "Original Spectrum"),
             col = c("green", "black", "gray","blue"), pt.cex=3, lty = c(1,4, 4,1),lwd = c(3,3,3,3), cex = 0.8)
      if ((sno == i || !(i%%n)) && saveFig) {
        if (file.exists(outpdf2) && !overwriteFig)
          cat("Can't save figure, file", outpdf2, "already exists.\n")
        else
          df = dev.copy2pdf(device=x11, file = outpdf2)
      }
    }
  } else {
    cat("No results found.\n")
  }
  warnRead<-options(warn = warnDef)
}
