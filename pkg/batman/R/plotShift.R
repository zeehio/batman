plotShift<-function(BM, metaName, plotHist = FALSE, breaks, perMult = FALSE, saveFig = TRUE, 
                    saveFigDir = BM$outputDir, prefixFig, overwriteFig = FALSE, showPlot)
{ 
  ## written by Dr. Jie Hao, Imperial College London
  ## plot multiplet shift posteriors in boxplot or histogram
  if (missing(BM))
    return(cat("Please input batman data list.\n"))
  warnDef<-options("warn")$warn
  warnRead<-options(warn = -1)
  ptype = "pdf"
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
  
  mName<-NULL  ## metabolites names
  mValue<-NULL   ## corresponding ppm values
  mLabel<-NULL    ## label to be shown in plot
  mind<-NULL  ## matching multiplets index
  m<-NULL
  m<-names(BM$deltaSam)
  if (!is.null(m)) 
  { 
    mN<- strsplit(m, "_")
  } else {
    return(cat("No metabolite to match.\n")) 
  }
  
  ns<-5     
  for (i in 1:length(mN))
  {
    mName<-rbind(mName, gsub("[:.:]","-",mN[[i]][1]))
    mValue<-rbind(mValue,mN[[i]][2])
  }  
  
  mLabel<-paste(mName,mValue,sep=" ")
  mind <-NULL
  
  if (!missing(metaName)) {  
    meta<-gsub(" ","-",metaName)
    mind<-which(!is.na(match(tolower(mName),tolower(meta))))    
  } else {
    metaName<-NULL
  }   
  
  fc2<-length(mind)
  if (!is.null(BM$deltaSam)) 
  {
    f<-nrow(BM$deltaSam)
    fc<-ncol(BM$deltaSam)
    
    sno<-BM$specRange
    ind<-fc/length(sno)   ## total number of multiplets
    ind2<-fc2/length(sno)  ## total number of matching multiplets
    outpdf1<-NULL
    perSpec <- NULL       
    plotAllOP<-NULL
    if (!perMult || length(sno) == 1)
    {
      for (i in 1:length(sno)) 
      {
        isno <- sno[i]
        if (!missing(prefixFig))
          outpdf1 <- paste(saveFigDir, "/", prefixFig,"_spec_",isno, "_",metaName,"_ppmShift.",ptype, sep="")
        else
          outpdf1 <- paste(saveFigDir,"/spec_",isno, "_",metaName,"_ppmShift.",ptype, sep="")
        #if (!is.null(metaName))
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
        if (missing(breaks))
        {
          breaks <- length(BM$deltaSam[,1])%/%3
          if (breaks <=0)
            breaks <- length(BM$deltaSam[,1])
        } 
        ## check if a metabolite name is matched, if not plotting all
        if (length(mind)>0) 
        {
          if (plotHist) 
          {
            ## plot result in histogram
            n<-length(mind[((i-1)*ind2+1):(i*ind2)])
            if (n <=5)
              par(mfrow=c(n, 1),oma = c(0, 0, 3, 0))
            else if (n <= 16 )
              par(mfrow=c(ceiling(n/4),4), oma = c(0, 0, 3, 0))
            else 
              par(mfrow=c(ceiling(n/5),5), oma = c(0, 0, 3, 0))
            for (j in 1:n)
            {
              hist(BM$deltaSam[,mind[((i-1)*ind2+1)+j-1]],col="gray",
                   main = paste(metaName, " at ", mValue[mind[((i-1)*ind2+1)+j-1]], " ppm",sep=""),
                   xlab = "ppm shift", breaks = breaks, prob=TRUE)
              lines(density(BM$deltaSam[,mind[((i-1)*ind2+1)+j-1]])) 
              v<-quantile(BM$deltaSam[,mind[((i-1)*ind2+1)+j-1]],p = c(2.5,97.5)/100, type = 3)
              abline(v = v, col = "red", lty=2, lwd = 2)
              legend("topright", legend = "2.5% & 97.5% quantiles", col = "red", pt.cex=2, lty =  2, cex = 0.8)
              box()
            }
            title(main = list(paste("NMR Spectrum ",isno,": ", BM$specTitle[2,i], sep=""),cex=1.5,font=3), outer=TRUE)                  
          } else {
            ## plot result in boxplot
            boxplot(as.data.frame(BM$deltaSam[,mind[((i-1)*ind2+1):(i*ind2)]]), 
                    main = paste("NMR Spectrum ", isno, ": ", BM$specTitle[2,i], "\n", metaName, sep=""), 
                    names = as.character(mValue[mind[((i-1)*ind2+1):(i*ind2)]]),
                    ylab="ppm shift",xaxt="n") 
            axis(1,at=1:length(mind[((i-1)*ind2+1):(i*ind2)]),adj=1,padj=0.5,
                 labels=as.character(mValue[mind[((i-1)*ind2+1):(i*ind2)]]),las=2) 
          }
        } else {
          if (plotHist) 
          {
            if (ind <=5)
              par(mfrow=c(ind, 1),oma = c(0, 0, 3, 0))
            else if (ind <= 16)
              par(mfrow=c(4, ceiling(ind/4)),oma = c(0, 0, 3, 0))
            else 
              par(mfrow=c(5, ceiling(ind/5)),oma = c(0, 0, 3, 0))
            for (j in 1:ind)
            {
              hist(BM$deltaSam[,((i-1)*ind+1)+j-1],col="gray",
                   main = paste(mName[((i-1)*ind+1)+j-1], " at ", mValue[((i-1)*ind+1)+j-1], " ppm",sep=""),
                   xlab = "ppm shift", breaks = breaks, prob=TRUE)
              lines(density(BM$deltaSam[,((i-1)*ind+1)+j-1]))
              v<-quantile(BM$deltaSam[,((i-1)*ind+1)+j-1],p = c(2.5,97.5)/100, type = 3)
              abline(v = v, col = "red", lty=2, lwd = 2)
              legend("topright", legend = "2.5% & 97.5% quantiles", col = "red", pt.cex=2, lty =  2, cex = 0.8)
              box()
            }
            title(main = list(paste("NMR Spectrum ",isno,": ", BM$specTitle[2,i], sep=""),cex=1.5,font=3), outer=TRUE)   
          } else {
            par(oma = c(4, 0, 0, 0))  
            boxplot(as.data.frame(BM$deltaSam[,((i-1)*ind+1):(i*ind)]), 
                    main = paste("NMR Spectrum ", isno, ": ", BM$specTitle[2,i], "\n", metaName, sep=""), 
                    names = as.character(mLabel[((i-1)*ind+1):(i*ind)]), ylab="ppm shift", xaxt="n")
            axis(1,at=1:length(mLabel[((i-1)*ind+1):(i*ind)]),adj=1,padj=0.5,
                 labels=as.character(mLabel[((i-1)*ind+1):(i*ind)]), las=2) 
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
          
          ##if (file.exists(outpdf1) && !overwriteFig)
          ##  cat("Can't save figure, file", outpdf1, "already exists.\n")
          ## else
          ##  df = dev.copy2pdf(device=x11, file = outpdf1)
        }
      }
    } else {
      ##  plot per multiplet
      if (length(mind)>0)
      {
        n<-length(mind[((1-1)*ind2+1):(1*ind2)])
        mLab<-mLabel[mind]
      } else {
        #cat("No matching metabolite found for spectrum, ploting all metabolites.\n")
        n<-ind
        mLab<-mLabel
      }
      for (j in 1:n) 
      {
        if (!missing(prefixFig))
          outpdf1 <- paste(saveFigDir, "/", prefixFig,"_",mLab[j],"_ppmShift.",ptype, sep="")
        else
          outpdf1 <- paste(saveFigDir,"/",mLab[j],"_ppmShift.",ptype, sep="") 
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
        
        if (missing(breaks))
        {
          breaks <- length(BM$deltaSam[,1])%/%3
          if (breaks <=0)
            breaks <- length(BM$deltaSam[,1])
        } 
        ## check if a metabolite name is matched, if not plotting all
        if (length(mind)>0) 
        {
          ## plot result in histogram
          if (plotHist) 
          {
            #n<-length(BM$deltaSam[,mind[((i-1)*ind2+1):(i*ind2)]])
            if (length(sno) <=5)
              par(mfrow=c(length(sno), 1),oma = c(0, 0, 3, 0))
            else if (length(sno) <= 16 )
              par(mfrow=c(ceiling(length(sno)/4),4), oma = c(0, 0, 3, 0))
            else 
              par(mfrow=c(ceiling(length(sno)/5),5), oma = c(0, 0, 3, 0))
            for (i in 1:length(sno))
            {
              isno <- sno[i]
              hist(BM$deltaSam[,mind[((i-1)*ind2+1)+j-1]],col="gray",
                   main = paste("NMR Spectrum ",isno,": ", BM$specTitle[2,i], sep=""),
                   xlab = "ppm shift", breaks = breaks,prob=TRUE)
              lines(density(BM$deltaSam[,mind[((i-1)*ind2+1)+j-1]]))             
              v<-quantile(BM$deltaSam[,mind[((i-1)*ind2+1)+j-1]],p = c(2.5,97.5)/100, type = 3)
              abline(v = v, col = "red", lty=2, lwd = 2)      
              legend("topright", legend = "2.5% & 97.5% quantiles",
                     col = "red", pt.cex=2, lty =  2, cex = 0.8)
              box()
            }
            title(main = list(paste(metaName, " at ", mValue[mind[((1-1)*ind2+1)+j-1]]," ppm",sep=""),cex=1.5,font=3), outer=TRUE)                  
          } else {
            mtmp<-BM$deltaSam[,mind[seq(1,length(mind),ind2)+j-1]]
            boxplot(as.data.frame(mtmp), 
                    main = paste(metaName, " at ", mValue[mind[((1-1)*ind2+1)+j-1]], " ppm",sep=""), 
                    names = as.character(t(BM$specTitle[2,])),
                    ylab="ppm shift", xaxt="n") 
            axis(1,at=1:sno,adj=1,padj=0.5,labels=as.character(t(BM$specTitle[2,])),las=2) 
          }
        } else {
          #cat("No matching metabolite found for spectrum", i, "ploting all multiplets.\n")
          if (plotHist) 
          {
            if (length(sno) <=5)
              par(mfrow=c(length(sno), 1),oma = c(0, 0, 3, 0))
            else if (length(sno) <= 16 )
              par(mfrow=c(ceiling(length(sno)/4),4),oma = c(0, 0, 3, 0))
            else 
              par(mfrow=c(ceiling(length(sno)/5),5),oma = c(0, 0, 3, 0))
            for (i in 1:length(sno))
            {
              isno <- sno[i]
              hist(BM$deltaSam[,((i-1)*ind+1)+j-1],col="gray",
                   main = paste("NMR Spectrum ",isno,": ", BM$specTitle[2,i], sep=""),
                   xlab = "ppm shift", breaks = breaks,prob=TRUE)
              lines(density(BM$deltaSam[,((i-1)*ind+1)+j-1]))                            
              v<-quantile(BM$deltaSam[,((i-1)*ind+1)+j-1],p = c(2.5,97.5)/100, type = 3)
              abline(v = v, col = "red", lty=2, lwd = 2)
              legend("topright", legend = "2.5% & 97.5% quantiles",
                     col = "red", pt.cex=2, lty =  2, cex = 0.8)
              box()
            }
            title(main = list(paste(mLab[j], " ppm",sep=""),cex=1.5,font=3), outer=TRUE)  
          } else {
            ## per Multiplet and ploting boxplot
            par(oma = c(4, 0, 0, 0))	                 
            boxplot(as.data.frame(BM$deltaSam[,(seq(1,fc,ind)+j-1)]), 
                    main = paste(mName[j], " at ", mValue[j], " ppm",sep=""), 
                    names = as.character(t(BM$specTitle[2,])), ylab="ppm shift", las=2,xaxt="n") 
            axis(1,at=1:length(sno),adj=1,padj=0.5,labels=as.character(t(BM$specTitle[2,])),las=2)                       
          }
        }
        #title(main = list(paste("ppm shift\nNMR Spectrum ",i,": ", BM$specTitle[2,i], sep=""),cex=1.5,font=3), outer=TRUE)                  
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
          
          #if (file.exists(outpdf1) && !overwriteFig)
          #  cat("Can't save figure, file", outpdf1, "already exists.\n")
          #else
          #  df = dev.copy2pdf(device=x11, file = outpdf1)
        }
      }
    }
  }
  warnRead<-options(warn = warnDef)
}