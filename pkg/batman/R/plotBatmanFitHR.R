plotBatmanFitHR<-function(BM, xfrom, xto, yfrom, yto, metaName, saveFig = TRUE, 
                          saveFigDir = BM$outputDir, prefixFig, rerun = FALSE, overwriteFig = FALSE)
{      
  ## written by Dr. Jie Hao, Imperial College London
  ## plot batman metabolite fitting results in its original resolution 
  if (missing(BM))
    return(cat("Please input batman data list.\n"))
  warnDef<-options("warn")$warn
  warnRead<-options(warn = -1)
  ## save in pdf format
  ptype = "pdf"
  cex = 0.8
  ns<-5
  nsH<-3
  ## set default values if missing input
  m<-row.names(BM$beta)
  if (!missing(metaName))
  {
    mind<-which(!is.na(match(tolower(m),tolower(metaName)))) 
    if (length(mind)==0) 
      cat("No matching metabolite found...\n")   
  } else {
    metaName <-NULL
    mind<-NULL
  }
  if (missing(xfrom))
    xfrom <- min(BM$sFitHR[1,1],BM$sFitHR[nrow(BM$sFitHR),1])
  if (missing(xto))
    xto <- max(BM$sFitHR[1,1],BM$sFitHR[nrow(BM$sFitHR),1])  
  
  if (xfrom > xto)
  {
    temp<-xfrom
    xfrom<-xto
    xto<-temp
  }
  
  if (missing(yfrom)) { yfrom <- 0 }
  if (missing(yto)) { ytoIP <- NULL} else {ytoIP <- yto}
  
  pindH<-which(BM$sFitHR[,1]<=xto & BM$sFitHR[,1]>=xfrom)
  if (length(pindH)== 0)
    pindH<-which(BM$sFitHR[,1]<=xfrom & BM$sFitHR[,1]>=xto)
  
  pind<-which(BM$sFit[,1]<=xto & BM$sFit[,1]>=xfrom)
  if (length(pind)== 0)
    pind<-which(BM$sFit[,1]<=xfrom & BM$sFit[,1]>=xto)
  
  sno<-BM$specRange
  n <- 2
  metaTmplty <-nsH
  metaTmplwd <-3
  plotcol<- sample(rainbow(nrow(BM$beta)))
  outpdf1<-NULL
  ## plot batman results
  if (!is.null(BM$sFitHR) && !rerun)
  {
    for (j in 1:length(sno))
    {
      jsno <- sno[j]
      if (!missing(prefixFig))
        outpdf1 <- paste(saveFigDir, "/", prefixFig,"_specFitHR_",jsno, "_",metaName,".",ptype, sep="")
      else
        outpdf1 <- paste(saveFigDir,"/specFitHR_",jsno, "_",metaName,".",ptype, sep="")
      x11(15,7)
      
      i = ns*(j-1)+1
      iH = nsH*(j-1)+1
      
      if (is.null(ytoIP)) { yto <- max(max(BM$sFitHR[pindH,iH+1]),BM$sFitHR[pindH,iH+2])}
      
      if (yfrom > yto)
      {
        temp<-yfrom
        yfrom<-yto
        yto<-temp
      }
      ## whole spectrum
      plot(BM$sFitHR[pindH,iH],BM$sFitHR[pindH,iH+1],type="l",xlim=rev(range(BM$sFitHR[pindH,iH])),xlab="ppm",
           ylab="Standardized Intensity", main=paste("NMR Spectrum ",jsno,": ",BM$specTitle[2,j], sep=" "), 
           ylim = c(yfrom, yto), lwd = 1, col = 4, lty = 1)
      points(BM$sFit[pind,i],BM$sFit[pind,i+1],col=3, lwd = 1)
      lines(BM$sFitHR[pindH,iH],BM$sFitHR[pindH,iH+2],col=2, lwd = 1, lty = 1)
      points(BM$sFit[pind,i],BM$sFit[pind,i+2],col=1, lwd = 1)
      ## add on named metabolite
      if (!missing(metaName) && length(mind)!= 0 )
      {
        lines(BM$sFitHR[pindH,iH],BM$beta[mind,j]*BM$metaTempHR[pindH,mind+(j-1)*nrow(BM$beta)],
              col=5, lwd = metaTmplwd, lty = metaTmplty )
        points(BM$sFit[pind,i],BM$beta[mind,j]*BM$metaTemp[pind,mind+(j-1)*nrow(BM$beta)],
               col=6, lwd = metaTmplwd)
        
        legend("topright", c("Original Spectrum (High Resolution)", "Original Spectrum (Downsampled)", "Metabolites Fit (High Resolution)",
                             "Metabolites Fit (Downsampled)", paste(m[mind]," (High Resolution)",sep = ""), paste(m[mind]," (Downsampled)",sep = "")), 
               pch = c(-1, 1, -1, 1, -1, 1), col=c(4,3,2,1,5,6), cex = cex,
               lty=c(1,-1,1,-1,metaTmplty, -1),lwd = c(rep(1,4),rep(metaTmplwd, 2)))
      }  else   {
        legend("topright", c("Original Spectrum (High Resolution)", "Original Spectrum (Downsampled)", 
                             "Metabolites Fit (High Resolution)", "Metabolites Fit (Downsampled)"), pch = c(-1,1,-1,1), col=c(4,3,2,1), 
               lty=c(1,-1,1,-1),lwd = c(rep(1,4)), cex = cex)
      }
      if (saveFig) {
        if (file.exists(outpdf1) && !overwriteFig)
          cat("Can't save figure, file", outpdf1, "already exists.\n")
        else
          df = dev.copy2pdf(device=x11, file = outpdf1)
      }
    }
  }
  ## plot rerun batman results 
  else if (!is.null(BM$sFitRerunHR) && rerun) 
  {     
    outpdf2 <-NULL
    for (j in 1:length(sno))
    {
      if (!missing(prefixFig))
        outpdf2 <- paste(saveFigDir, "/", prefixFig,"_specFitRerunHR_",jsno, "_",metaName,".",ptype, sep="")
      else
        outpdf2 <- paste(saveFigDir,"/specFitRerunHR_",jsno, "_",metaName,".",ptype, sep="")
      x11(15,7)
      
      i = ns*(j-1)+1
      iH = nsH*(j-1)+1
      
      if (is.null(ytoIP))
        yto <- max(max(BM$sFitRerunHR[pindH,iH+1]),BM$sFitRerunHR[pindH,iH+2])
      
      if (yfrom > yto)
      {
        temp<-yfrom
        yfrom<-yto
        yto<-temp
      }
      ## whole spectrum
      plot(BM$sFitRerunHR[pindH,iH],BM$sFitRerunHR[pindH,iH+1],type="l",xlim=rev(range(BM$sFitRerunHR[pindH,iH])),xlab="ppm",
           ylab="Standardized Intensity", main=paste("NMR Spectrum ",jsno, ": ", BM$specTitle[2,j], " (Rerun)", sep=""),
           ylim = c(yfrom, yto), lwd = 1, col = 4, lty = 1)
      points(BM$sFitRerun[pind,i],BM$sFitRerun[pind,i+1],col=3, lwd = 1, lty = 1, type = "o")
      lines(BM$sFitRerunHR[pindH,iH],BM$sFitRerunHR[pindH,iH+2],col=2, lwd = 1, lty = 1)
      points(BM$sFitRerun[pind,i],BM$sFitRerun[pind,i+2],col=1, lwd = 1, lty = 1, type = "o")
      # add on named metabolite 
      if (!missing(metaName) && length(mind)!= 0 )
      {
        lines(BM$sFitRerunHR[pindH,iH],BM$betaRerun[mind,j]*BM$metaTempRerunHR[pindH,mind+(j-1)*nrow(BM$betaRerun)],
              col=5, lwd = metaTmplwd, lty = metaTmplty )
        lines(BM$sFitRerun[pind,i],BM$betaRerun[mind,j]*BM$metaTempRerun[pind,mind+(j-1)*nrow(BM$betaRerun)],
              col=6, lwd = metaTmplwd, lty = metaTmplty, type = "o")
        
        legend("topright", c("Original Spectrum (High Resolution)", "Original Spectrum (Downsampled)",
                             "Metabolites Fit (High Resolution)", "Metabolites Fit (Downsampled)",
                             paste(m[mind]," (High Resolution)",sep = ""), paste(m[mind]," (Downsampled)",sep = "")), 
               pch = c(-1, 1, -1, 1, -1, 1), col=c(4,3,2,1,5,6), ncol = 2, cex = cex,
               lty=c(1,-1,1,-1,metaTmplty, -1),lwd = c(rep(1,4),rep(metaTmplwd, 2)))
      } else {
        legend("topright", c("Original Spectrum (High Resolution)", "Original Spectrum (Downsampled)", 
                             "Metabolites Fit (High Resolution)", "Metabolites Fit (Downsampled)"), pch = c(-1,1,-1,1), col=c(4,3,2,1), 
               lty=c(1,-1,1,-1),lwd = c(rep(1,4)), cex = cex)
      }
      if (saveFig) {
        if (file.exists(outpdf2) && !overwriteFig)
          cat("Can't save figure, file", outpdf2, "already exists.\n")
        else
          df = dev.copy2pdf(device=x11, file = outpdf2)
      }
    }
  } else {
    cat("No high resolution results found.\n")
  }
  warnRead<-options(warn = warnDef)
}
