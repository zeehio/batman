plotChemShiftDist<-function(BM, metaName, breaks = 20, xlim, 
                            saveFig = TRUE, saveFigDir = BM$outputDir,
                            prefixFig, overwriteFig = FALSE, showPlot = TRUE)
{
  ## written by Dr. Jie Hao, Imperial College London
  ## Histogram of the mean posterior estimated chemical shifts for the multiplets of  
  ## a given metabolite across a series of spectra. 
  
  setxlim <- FALSE
  pdfdev = FALSE
  
  if (!missing(xlim))
  {
    setxlim <- TRUE
  }
  
  if (missing(metaName))
  {
  metaName <- rownames(BM$beta)
  }
  multiName <- rownames(BM$delta)
  
  myVector <- strsplit(multiName, "_")
  multiName1 <-NULL
  for (i in 1:length(myVector))
  {
    multiName1[i] <- myVector[[i]][1]
  }
  
  multiName1 <- gsub("[:.:]", " ", multiName1)
  metaName <- gsub("[:.:]", " ", metaName)
  metaName <- gsub("-", " ", metaName)
  
  multiNameLeg <- colnames(BM$deltaSam)[1:length(multiName1)]
  multiNameLeg <- gsub("_", " ", multiNameLeg)
  nMeta <- length(metaName)
  
  for (i2 in 1:nMeta)
  {
    outpdf1 <- paste(saveFigDir, "/chemShiftDist_", metaName[i2], ".pdf", sep="")  
    if ((!showPlot && overwriteFig) || (!showPlot && (!file.exists(outpdf1))))
    {
      pdf(outpdf1)  
      pdfdev = TRUE
    }        
    else if (!showPlot && (file.exists(outpdf1) && !overwriteFig))
      cat("Can't save figure, file", outpdf1, "already exists.\n")
    else
      x11()
    
    mid <- which(metaName[i2] == multiName1) 
    if (length(mid)>25) {
      cat("How many multiplets does", metaName[i2], "have? (Haven't program it to show more than 25)\n" )
    } else if (length(mid) > 9) {
      par(mfrow=c(5,5), oma = c(0, 0, 3, 0))    
    } else if (length(mid)>4) {
      par(mfrow=c(3,3), oma = c(0, 0, 3, 0))
    } else if (length(mid) >1){
      par(mfrow=c(2,2), oma = c(0, 0, 3, 0))
    }
    for (i in mid)
    { 
      if (setxlim)
      {
        hist(t(BM$delta[i,]),col="gray", main = multiNameLeg[i], xlim,
             xlab = expression(paste(Delta,"chemical shift")), breaks)  
      } else {
        hist(t(BM$delta[i,]),col="gray", main = multiNameLeg[i], 
             xlab = expression(paste(Delta,"chemical shift")), breaks)  
      }
    }
    if (saveFig) 
    {
      if (pdfdev)
      {
        pdfoff = dev.off()    
        pdfdev = FALSE
      }
      else if (showPlot && (file.exists(outpdf1) && !overwriteFig))
        cat("Can't save figure, file", outpdf1, "already exists.\n")
      else
        df = dev.copy2pdf(device=x11, file = outpdf1)
      
      ## if (file.exists(outpdf1) && !overwriteFig)
      ##   cat("Can't save figure, file", outpdf1, "already exists.\n")
      ## else
      ##   df = dev.copy2pdf(device=x11, file = outpdf1)
    }
  }
}
