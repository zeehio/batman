plotChemShiftDist<-function(BM, metaName, breaks = 20, xlim, 
                            saveFig = TRUE, saveFigDir = BM$outputDir,
                            prefixFig, overwriteFig = FALSE, showPlot)
{
  ## written by Dr. Jie Hao, Imperial College London
  ## Histogram of the mean posterior estimated chemical shifts for the multiplets of  
  ## a given metabolite across a series of spectra. 
  
  setxlim <- FALSE
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
    {  
      cat("\nCan't save figure, file", outpdf1, "already exists.\n")
      tmpOP <- strsplit(outpdf1, "[.]")
      outpdf1 <- paste(tmpOP[[1]][1], "_", format(Sys.time(), "%d_%b_%H_%M_%S"), ".", tmpOP[[1]][2], sep = "")
      cat("Figure saved in new file \"", outpdf1, "\".\n")
      pdf(outpdf1)  
      pdfdev = TRUE
    } 
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
      ##   df = dev.copy2pdf(device=x11, file = outpdf1)
    }
  }
}
