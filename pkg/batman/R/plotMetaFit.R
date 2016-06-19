plotMetaFit<-function(BM, from, to, metaName, saveFig = TRUE,saveFigDir = BM$outputDir, 
                      prefixFig, rerun = FALSE, overwriteFig = FALSE, showPlot)
{
  ## written by Dr. Jie Hao, Imperial College London
  ## plot metabolites fit posteriors with 95% credible Interval
  if (missing(BM))
    return(cat("Please input batman data list.\n"))
  
  warnDef<-options("warn")$warn
  warnRead<-options(warn = -1)
  ptype = "pdf"
  n <- 2
  ns<-5
  pdfdev = FALSE
  
  ## os information
  os <- NULL
  if (missing(showPlot))
  {
    sysinf <- Sys.info()
    os <- "notlisted"
    if (!is.null(sysinf)){
      os1 <- sysinf['sysname']
      if (os1 == 'Darwin')
      {os <- "osx"}
      else if (grepl("windows", tolower(os1)))
      {os<- "win"}       
      else if (grepl("linux", tolower(os1)))
      {os<- "linux" }
    } else { ## mystery machine
      #os <- .Platform$OS.type
      if (grepl("^darwin", R.version$os))
        os <- "osx"
      if (grepl("linux-gnu", R.version$os))
        os <- "linux"
    }
  }
  
  if (!is.null(os))
  {
    if (os == 'win' || os == 'osx')
    { showPlot <- TRUE }
    else 
    { #if (os == 'linux')
      showPlot <- FALSE
      cat("\nThis operating system may not support X11, no plot will be displayed, figures in .pdf format will be saved in output folder.")
      cat("\nCheck input argument 'showPlot' for more detail.")
    }
  }
  
  
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
          
          if ((!showPlot && overwriteFig) || (!showPlot && (!file.exists(outpdf1))))
          {
              pdf(outpdf1,15,7)  
              pdfdev = TRUE
          }         
          else if (!showPlot && (file.exists(outpdf1) && !overwriteFig))
          {  
            cat("\nCan't save figure, file", outpdf1, "already exists.\n")
            tmpOP <- strsplit(outpdf1, "[.]")
            outpdf1 <- paste(tmpOP[[1]][1], "_", format(Sys.time(), "%d_%b_%H_%M_%S"), ".", tmpOP[[1]][2], sep = "")
            cat("Figure saved in new file \"", outpdf1, "\".\n")
            pdf(outpdf1,15,7)  
            pdfdev = TRUE
          }
          else
            x11(15,7)
          par(mfrow=c(n,1))  
        } else {
          if (!missing(prefixFig))
            outpdf1 <- paste(saveFigDir, "/", prefixFig, "_spec_",sno[i],"_mFitSam_", metaName, ".",ptype, sep="")
          else
            outpdf1 <- paste(saveFigDir,"/spec_",sno[i],"_mFitSam_", metaName, ".",ptype, sep="")
          
          if ((!showPlot && overwriteFig) || (!showPlot && (!file.exists(outpdf1))))
          {
              pdf(outpdf1,15,7)  
              pdfdev = TRUE
          }           
          else if (!showPlot && (file.exists(outpdf1) && !overwriteFig))
          {  
            cat("\nCan't save figure, file", outpdf1, "already exists.\n")
            tmpOP <- strsplit(outpdf1, "[.]")
            outpdf1 <- paste(tmpOP[[1]][1], "_", format(Sys.time(), "%d_%b_%H_%M_%S"), ".", tmpOP[[1]][2], sep = "")
            cat("Figure saved in new file \"", outpdf1, "\".\n")
            pdf(outpdf1,15,7)  
            pdfdev = TRUE
          }
          else
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
        v<-apply(metaFittemp, 1, quantile, probs = c(2.5,97.5)/100,  na.rm = TRUE, type = 3) 
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
        if (pdfdev)
        {
          pdfoff = dev.off()    
          pdfdev = FALSE
        }
        else if (showPlot && (file.exists(outpdf1) && !overwriteFig))
        {  
          cat("\nCan't save figure, file", outpdf1, "already exists.\n")
          tmpOP <- strsplit(outpdf1, "[.]")
          outpdf1 <- paste(tmpOP[[1]][1], "_", format(Sys.time(), "%d_%b_%H_%M_%S"), ".", tmpOP[[1]][2], sep = "")
          cat("Figure saved in new file \"", outpdf1, "\".\n")
          df = dev.copy2pdf(device=x11, file = outpdf1)
        }
        else
          df = dev.copy2pdf(device=x11, file = outpdf1)
        
        ## if (file.exists(outpdf1) && !overwriteFig)
        ##   cat("Can't save figure, file", outpdf1, "already exists.\n")
        ## else
        ##    df = dev.copy2pdf(device=x11, file = outpdf1)
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
          
          if ((!showPlot && overwriteFig) || (!showPlot && (!file.exists(outpdf2))))
          {
              pdf(outpdf2,15,7)  
              pdfdev = TRUE
          }            
          else if (!showPlot && (file.exists(outpdf2) && !overwriteFig))
          {  
            cat("\nCan't save figure, file", outpdf2, "already exists.\n")
            tmpOP <- strsplit(outpdf2, "[.]")
            outpdf2 <- paste(tmpOP[[1]][1], "_", format(Sys.time(), "%d_%b_%H_%M_%S"), ".", tmpOP[[1]][2], sep = "")
            cat("Figure saved in new file \"", outpdf2, "\".\n")
            pdf(outpdf2,15,7)  
            pdfdev = TRUE
          } 
          else
            x11(15,7)
          par(mfrow=c(n,1))	
        } else {
          if (!missing(prefixFig))
            outpdf2 <- paste(saveFigDir, "/", prefixFig, "_specRerun_",sno[i],"_mFitSam_", metaName,".",ptype, sep="")
          else
            outpdf2 <- paste(saveFigDir,"/specRerun_",sno[i],"_mFitSam_", metaName,".",ptype, sep="")
          
          if ((!showPlot && overwriteFig) || (!showPlot && (!file.exists(outpdf2))))
          {
              pdf(outpdf2,15,7)  
              pdfdev = TRUE
          }          
          else if (!showPlot && (file.exists(outpdf2) && !overwriteFig))
          {  
            cat("\nCan't save figure, file", outpdf2, "already exists.\n")
            tmpOP <- strsplit(outpdf2, "[.]")
            outpdf2 <- paste(tmpOP[[1]][1], "_", format(Sys.time(), "%d_%b_%H_%M_%S"), ".", tmpOP[[1]][2], sep = "")
            cat("Figure saved in new file \"", outpdf2, "\".\n")
            pdf(outpdf2,15,7)  
            pdfdev = TRUE
          } 
          else
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
        v<-apply(metaFittemp, 1, quantile, probs = c(2.5,97.5)/100,  na.rm = TRUE, type = 3) 
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
        if (pdfdev)
        {
          pdfoff = dev.off()    
          pdfdev = FALSE
        }       
        else if (showPlot && (file.exists(outpdf2) && !overwriteFig))
        {  
          cat("\nCan't save figure, file", outpdf2, "already exists.\n")
          tmpOP <- strsplit(outpdf2, "[.]")
          outpdf2 <- paste(tmpOP[[1]][1], "_", format(Sys.time(), "%d_%b_%H_%M_%S"), ".", tmpOP[[1]][2], sep = "")
          cat("Figure saved in new file \"", outpdf2, "\".\n")
          df = dev.copy2pdf(device=x11, file = outpdf2)
        }
        else
          df = dev.copy2pdf(device=x11, file = outpdf2)
        
        #if (file.exists(outpdf2) && !overwriteFig)
        #  cat("Can't save figure, file", outpdf2, "already exists.\n")
        # else
        #  df = dev.copy2pdf(device=x11, file = outpdf2)
      }
    }
  } else {
    cat("No results found.\n")
  }
  warnRead<-options(warn = warnDef)
}
