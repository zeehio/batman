batmanrerun<- function(BM, figBatmanFit = TRUE, listMeta = FALSE, 
                       figRelCon = FALSE, figMetaFit = FALSE)
{   
    ## rerun batman() with fixed ppm position for all multiplets
    ## directories for input parameters and data files
    dir1 <- paste(BM$outputDir,"/batmanOptions.txt",sep="")
   	dir2<-paste(BM$outputDir,"/NMRdata.txt",sep="")
	dir3<-paste(BM$outputDir,"/metabolitesList.txt",sep="")
	dir4<-paste(BM$outputDir,"/multi_data.dat",sep="")
    
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
	filedir<-c(dir1,dir2,dir3,dir4,dir5)
	
	## read in parameters from batmanOption.txt 
    con  <- file(dir1, open = "r")
    oneLine <- readLines(con, n = 30, warn = FALSE)
    fL<-substr(oneLine,1,1)
    nL<-which(is.na(match(fL,"%")))
    myVector <- strsplit(oneLine[nL[2]], ":")
    sno <- as.numeric(myVector[[1]][2])
    myVector <- strsplit(oneLine[nL[8]], ":")
    itoBI <- as.numeric(myVector[[1]][2])
    myVector <- strsplit(oneLine[nL[10]], ":")
    fixeff <- as.numeric(myVector[[1]][2])
    myVector <- strsplit(oneLine[nL[11]], ":")
    itoRr <- as.numeric(myVector[[1]][2])
    close(con)
    
    cat("Rerunning batman for ", itoRr, " iterations.\n")
    ## choose whether to parallelize between spectra
    if (sno>1 && fixeff == 0)    
    {
		cat("\nHow many parallel processes (multicores) do you want to run the multi-spectra analysis?")
		cat("\n(Enter 1 for running them sequentially.)\n")
		cat("\n Parallel processing of multi spectra currently cannot display the progress\n")
        cat(" bar (or any words), if you input is > 1, please be patient for the results :)\n\n")
		wr<- getinput(lowlim=1,highlim=20)  
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
                multispec1<- foreach (nospec = 1:sno,.packages="batman") %dopar% {
                    pBar <- txtProgressBar(min =0, max = itoRr, style = 3)
                    out<-.Call("batman",filedir,as.integer(bn),as.integer(itoRr),as.integer(rr),as.integer(nospec-1),pBar,PACKAGE = "batman")
                    close( pBar )
        }})[3]} else if (wr == 1 && fixeff == 1) {
        	stime1 <- system.time ({
            	multispec1<- for(nospec in 1:1) {
                    pBar <- txtProgressBar(min =0, max = itoRr, style = 3)
                	out<-.Call("batman",filedir,as.integer(bn),as.integer(itoRr),as.integer(rr),as.integer(sno-1),pBar,PACKAGE = "batman")
                	close( pBar ) 
        }})[3]}	else {
        	stime1 <- system.time ({
            	multispec1<- for(nospec in 1:sno) {
                    pBar <- txtProgressBar(min =0, max = itoRr, style = 3)
                	out<-.Call("batman",filedir,as.integer(bn),as.integer(itoRr),as.integer(rr),as.integer(nospec-1),pBar,PACKAGE = "batman")
                	close( pBar ) 
    	}})[3]}	
        cat (" For batman rerun, time ")
        print(stime1)
        cat ("seconds.\n")
        if (wr>1) {stopCluster(cl)}
    }

    cat("Reading in saved data in folder\n")
    ## read in results into R
	BMR<-readBatmanOutput(BM$outputDir) 
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

plotBatmanFitHR<-function(BM, xfrom, xto, yfrom, yto, metaName, saveFig = TRUE, 
saveFigDir = BM$outputDir, prefixFig, rerun = FALSE)
{      
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
    if (missing(yto)) { ytoIP <- NULL}
        
    pindH<-which(BM$sFitHR[,1]<=xto & BM$sFitHR[,1]>=xfrom)
    if (length(pindH)== 0)
        pindH<-which(BM$sFitHR[,1]<=xfrom & BM$sFitHR[,1]>=xto)

    pind<-which(BM$sFit[,1]<=xto & BM$sFit[,1]>=xfrom)
    if (length(pind)== 0)
        pind<-which(BM$sFit[,1]<=xfrom & BM$sFit[,1]>=xto)
   
    sno<-length(BM$sFitHR)/nsH
	n <- 2
    metaTmplty <-nsH
    metaTmplwd <-3
    plotcol<- sample(rainbow(nrow(BM$beta)))
    outpdf1<-NULL
    ## plot batman results
    if (!is.null(BM$sFitHR) && !rerun)
    {
	    for (j in 1:sno)
        {
            if (!missing(prefixFig))
                outpdf1 <- paste(saveFigDir, "/", prefixFig,"_specFitHR_",j, "_",metaName,".",ptype, sep="")
            else
                outpdf1 <- paste(saveFigDir,"/specFitHR_",j, "_",metaName,".",ptype, sep="")
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
        		 ylab="Standardized Intensity", main=paste("NMR Spectrum ",j,": ",BM$specTitle[2,j], sep=" "), 
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
            }	else   {
                legend("topright", c("Original Spectrum (High Resolution)", "Original Spectrum (Downsampled)", 
        		"Metabolites Fit (High Resolution)", "Metabolites Fit (Downsampled)"), pch = c(-1,1,-1,1), col=c(4,3,2,1), 
                lty=c(1,-1,1,-1),lwd = c(rep(1,4)), cex = cex)
            }
            if (saveFig) {
                if (file.exists(outpdf1))
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
		for (j in 1:sno)
        {
            if (!missing(prefixFig))
                outpdf2 <- paste(saveFigDir, "/", prefixFig,"_specFitRerunHR_",j, "_",metaName,".",ptype, sep="")
            else
                outpdf2 <- paste(saveFigDir,"/specFitRerunHR_",j, "_",metaName,".",ptype, sep="")
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
    			 ylab="Standardized Intensity", main=paste("NMR Spectrum ",j, ": ", BM$specTitle[2,j], " (Rerun)", sep=""),
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
                if (file.exists(outpdf2))
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

plotShift<-function(BM, metaName, plotHist = FALSE, breaks, perMult = FALSE, saveFig = TRUE, saveFigDir = BM$outputDir, prefixFig)
{ 
## plot multiplet shift posteriors in boxplot or histogram
    if (missing(BM))
	return(cat("Please input batman data list.\n"))
    warnDef<-options("warn")$warn
    warnRead<-options(warn = -1)
    ptype = "pdf"
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
		
        sno<-length(BM$sFit)/ns
        ind<-fc/sno   ## total number of multiplets
        ind2<-fc2/sno  ## total number of matching multiplets
        outpdf1<-NULL
        perSpec <- NULL       
        plotAllOP<-NULL
        if (!perMult || sno == 1)
        {
            for (i in 1:sno) 
            {
                if (!missing(prefixFig))
				outpdf1 <- paste(saveFigDir, "/", prefixFig,"_spec_",i, "_",metaName,"_ppmShift.",ptype, sep="")
                else
				outpdf1 <- paste(saveFigDir,"/spec_",i, "_",metaName,"_ppmShift.",ptype, sep="")
#if (!is.null(metaName))
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
						par(mfrow=c(ceiling(n/4)),4, oma = c(0, 0, 3, 0))
						else 
						par(mfrow=c(ceiling(n/5)),5, oma = c(0, 0, 3, 0))
						for (j in 1:n)
						{
							hist(BM$deltaSam[,mind[((i-1)*ind2+1)+j-1]],col="gray",
								 main = paste(metaName, " at ", mValue[mind[((i-1)*ind2+1)+j-1]], " ppm",sep=""),
								 xlab = "ppm shift", breaks = breaks, prob=TRUE)
							lines(density(BM$deltaSam[,mind[((i-1)*ind2+1)+j-1]])) 
							v<-quantile(BM$deltaSam[,mind[((i-1)*ind2+1)+j-1]],p = c(2.5,97.5)/100)
							abline(v = v, col = "red", lty=2, lwd = 2)
							legend("topright", legend = "2.5% & 97.5% quantiles", col = "red", pt.cex=2, lty =  2, cex = 0.8)
							box()
						}
						title(main = list(paste("NMR Spectrum ",i,": ", BM$specTitle[2,i], sep=""),cex=1.5,font=3), outer=TRUE)                  
					} else {
## plot result in boxplot
						boxplot(as.data.frame(BM$deltaSam[,mind[((i-1)*ind2+1):(i*ind2)]]), 
								main = paste("NMR Spectrum ", i, ": ", BM$specTitle[2,i], "\n", metaName, sep=""), 
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
							v<-quantile(BM$deltaSam[,((i-1)*ind+1)+j-1],p = c(2.5,97.5)/100)
							abline(v = v, col = "red", lty=2, lwd = 2)
							legend("topright", legend = "2.5% & 97.5% quantiles", col = "red", pt.cex=2, lty =  2, cex = 0.8)
							box()
						}
						title(main = list(paste("NMR Spectrum ",i,": ", BM$specTitle[2,i], sep=""),cex=1.5,font=3), outer=TRUE)   
					} else {
						par(oma = c(4, 0, 0, 0))	
						boxplot(as.data.frame(BM$deltaSam[,((i-1)*ind+1):(i*ind)]), 
								main = paste("NMR Spectrum ", i, ": ", BM$specTitle[2,i], "\n", metaName, sep=""), 
								names = as.character(mLabel[((i-1)*ind+1):(i*ind)]), ylab="ppm shift", xaxt="n")
						axis(1,at=1:length(mLabel[((i-1)*ind+1):(i*ind)]),adj=1,padj=0.5,
							 labels=as.character(mLabel[((i-1)*ind+1):(i*ind)]), las=2) 
					}
                }
				if (saveFig) 
				{
					if (file.exists(outpdf1))
					cat("Can't save figure, file", outpdf1, "already exists.\n")
					else
					df = dev.copy2pdf(device=x11, file = outpdf1)
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
						if (sno <=5)
						par(mfrow=c(sno, 1),oma = c(0, 0, 3, 0))
						else if (sno <= 16 )
						par(mfrow=c(ceiling(sno/4)),4, oma = c(0, 0, 3, 0))
						else 
						par(mfrow=c(ceiling(sno/5)),5, oma = c(0, 0, 3, 0))
						for (i in 1:sno)
						{
							hist(BM$deltaSam[,mind[((i-1)*ind2+1)+j-1]],col="gray",
								 main = paste("NMR Spectrum ",i,": ", BM$specTitle[2,i], sep=""),
								 xlab = "ppm shift", breaks = breaks,prob=TRUE)
							lines(density(BM$deltaSam[,mind[((i-1)*ind2+1)+j-1]]))             
							v<-quantile(BM$deltaSam[,mind[((i-1)*ind2+1)+j-1]],p = c(2.5,97.5)/100)
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
						if (sno <=5)
						par(mfrow=c(sno, 1),oma = c(0, 0, 3, 0))
						else if (sno <= 16 )
						par(mfrow=c(ceiling(sno/4)),4,oma = c(0, 0, 3, 0))
						else 
						par(mfrow=c(ceiling(sno/5)),5,oma = c(0, 0, 3, 0))
						for (i in 1:sno)
						{
							hist(BM$deltaSam[,((i-1)*ind+1)+j-1],col="gray",
								 main = paste("NMR Spectrum ",i,": ", BM$specTitle[2,i], sep=""),
								 xlab = "ppm shift", breaks = breaks,prob=TRUE)
							lines(density(BM$deltaSam[,((i-1)*ind+1)+j-1]))                            
							v<-quantile(BM$deltaSam[,((i-1)*ind+1)+j-1],p = c(2.5,97.5)/100)
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
						axis(1,at=1:sno,adj=1,padj=0.5,labels=as.character(t(BM$specTitle[2,])),las=2)                       
					}
				}
#title(main = list(paste("ppm shift\nNMR Spectrum ",i,": ", BM$specTitle[2,i], sep=""),cex=1.5,font=3), outer=TRUE)                  
				if (saveFig) 
				{
					if (file.exists(outpdf1))
					cat("Can't save figure, file", outpdf1, "already exists.\n")
                   else
                       df = dev.copy2pdf(device=x11, file = outpdf1)
               }
           }
       }
   }
   warnRead<-options(warn = warnDef)
}


plotMetaFit<-function(BM, from, to, metaName, saveFig = TRUE,saveFigDir = BM$outputDir, prefixFig, rerun = FALSE)
{
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
    
        sno<-length(BM$sFit)/ns
        ind<-fc/sno
        outpdf1<-NULL
        for (i in 1:sno) 
        {
            ## set subplot
            if ((i%%n) == 1){ 
                if ((sno-i)>=(n-1)) {
                    if (!missing(prefixFig))
                        outpdf1 <- paste(saveFigDir, "/", prefixFig,"_spec_",i,"to",i+n-1,"_mFitSam_", metaName, ".",ptype, sep="")
                    else
                        outpdf1 <- paste(saveFigDir,"/spec_",i,"to",i+n-1,"_mFitSam_", metaName, ".",ptype, sep="")
                    x11(15,7)
                    par(mfrow=c(n,1))	
                } else {
                    if (!missing(prefixFig))
                        outpdf1 <- paste(saveFigDir, "/", prefixFig, "_spec_",i,"_mFitSam_", metaName, ".",ptype, sep="")
                    else
                        outpdf1 <- paste(saveFigDir,"/spec_",i,"_mFitSam_", metaName, ".",ptype, sep="")
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
               maintitle <- paste("NMR Spectrum ", i, ": ",BM$specTitle[2,i], "\n", metaName, " fit sample", sep="")  
               metaFitFinal <-(BM$beta[mind,i]*BM$metaTemp[pind,mind+(i-1)*nrow(BM$beta)])
            } else {
               metaFittemp<-BM$metaFitSam[,((i-1)*ind+1):(i*ind)]
               maintitle <- paste("NMR Spectrum ", i, ": ",BM$specTitle[2,i], "\nMetabolite fit sample", sep="")
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
            if ((sno == i || !(i%%n)) && saveFig) {
                if (file.exists(outpdf1))
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
    
        sno<-length(BM$sFitRerun)/ns
        ind<-fc/sno
        outpdf2<-NULL
        for (i in 1:sno) {
            ## set subplot
            if ((i%%n) == 1){ 
                if ((sno-i)>=(n-1)) {
                    if (!missing(prefixFig))
                        outpdf2 <- paste(saveFigDir, "/", prefixFig,"_specRerun_",i,"to",i+n-1,"_mFitSam_", metaName, ".",ptype, sep="")
                    else
                        outpdf2 <- paste(saveFigDir,"/specRerun_",i,"to",i+n-1,"_mFitSam_", metaName, ".",ptype, sep="")
                    x11(15,7)
                    par(mfrow=c(n,1))	
                } else {
                    if (!missing(prefixFig))
                        outpdf2 <- paste(saveFigDir, "/", prefixFig, "_specRerun_",i,"_mFitSam_", metaName,".",ptype, sep="")
                    else
                        outpdf2 <- paste(saveFigDir,"/specRerun_",i,"_mFitSam_", metaName,".",ptype, sep="")
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
               maintitle <- paste("NMR Spectrum ", i, ": ", BM$specTitle[2,i], "(Rerun)\n", metaName, " fit sample", sep="")  
               metaFitFinal <-(BM$betaRerun[mind,i]*BM$metaTempRerun[pind,mind+(i-1)*nrow(BM$betaRerun)])
            } else {
                metaFittemp<-BM$metaFitSamRerun[,((i-1)*ind+1):(i*ind)]
                maintitle <- paste("NMR Spectrum ", i, ": ", BM$specTitle[2,i], "(Rerun)\nMetabolite fit sample", sep="")
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
                if (file.exists(outpdf2))
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

readBruker<-function(BrukerDataDir)
{
    warnDef<-options("warn")$warn
    warnRead<-options(warn = -1)
    datapath<-BrukerDataDir
    ## read in bruker spectra
    ## find the data files
	pfile <-list.files(path = datapath, pattern = "^procs$", all.files = FALSE,full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
	rfile <-list.files(path = datapath, pattern = "^1r$", all.files = FALSE,full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
	L<-length(pfile)
	Lr<-length(rfile)
	sa <- NULL
	snam <- NULL
	if (L==0 || Lr==0 || L!=Lr)
	{
    	return (cat("Bruker file does not exist in datapath, or other problems with bruker files...\n"))
	} else {
    	for (i in 1:L)
    	{    
            con  <- file(pfile[i], open = "r")
            aLine <- readLines(con, n = -1, warn = FALSE)
            myV <- strsplit(aLine, "=")    
            close(con)
        
            for (j in 1:length(myV))
            {
                if (match("##$OFFSET",myV[[j]][1],nomatch = 0))
                {    offset <- as.numeric(myV[[j]][2]);
                }
                if (match("##$SW_p",myV[[j]][1],nomatch = 0))
                {    sw <- as.numeric(myV[[j]][2]);
                }
                if (match("##$SF",myV[[j]][1],nomatch = 0))
                {
                    sf <- as.numeric(myV[[j]][2]);
                }
                if (match("##$SI",myV[[j]][1],nomatch = 0))
                {  
                    si <- as.numeric(myV[[j]][2]);
                }
                if (match("##$BYTORDP",myV[[j]][1],nomatch = 0))
                {    bytordp <- as.numeric(myV[[j]][2]);
                }
                if (match("##$NC_proc",myV[[j]][1],nomatch = 0))
                {
                    ncproc <- as.numeric(myV[[j]][2]);
                }
            }
        
            if (bytordp==0){machine_format =  "little"}
            else {machine_format = "big"}
        
            s<-readBin(rfile[i], what="int",70000, size = 4, signed = T, endian =machine_format)
            s<- ((2^ncproc)* s)
            nspec <- length(s)
        
            swp <- sw/sf
            dppm <- swp/(nspec-1)
            ppm<-offset
            ppm<-seq(offset,(offset-swp),by=-dppm)
            sa<- cbind(sa,s)
            ## find corresponding title
            stitle<-paste(substr(rfile[i],1,nchar(rfile[i])-2),"title",sep="")
            if (!file.exists(stitle))
            stitle<-paste(substr(rfile[i],1,nchar(rfile[i])-2),"TITLE",sep="")
            if (file.exists(stitle))
            {
                if (!file.info(stitle)$size == 0)
                {
                    con<-file(stitle,open="r")
                    ntem <- readLines(con, n = 1, warn = FALSE)
                    close(con)
                } else {
                    sT <- strsplit(rfile[i], "/")
                    sTitle <-sT[[1]]         
                    lsT<-length(sTitle)
                    if (lsT>4)
                    ntem<-paste(sTitle[lsT-4],"_",sTitle[lsT-3],"_",sTitle[lsT-1],sep="")
                    else if (lsT>3)
                    ntem<-paste(sTitle[lsT-3],"_",sTitle[lsT-1],sep="")
                    else if (lsT>=1)
                    ntem<-paste(sTitle[lsT-1],sep="")
                    else
                    ntem<-i
                }
            } else {
                sT <- strsplit(rfile[i], "/")
                sTitle <-sT[[1]]         
                lsT<-length(sTitle)
                if (lsT>4)
                ntem<-paste(sTitle[lsT-4],"_",sTitle[lsT-3],"_",sTitle[lsT-1],sep="")
                else if (lsT>3)
                ntem<-paste(sTitle[lsT-3],"_",sTitle[lsT-1],sep="")
                else if (lsT>=1)
                ntem<-paste(sTitle[lsT-1],sep="")
                else
                ntem<-i
            }
            snam<- cbind(snam, ntem)            
        }
    }
    snam <- cbind("ppm", snam)
    sa <- cbind(ppm,sa)
    colnames(sa)<- snam
    warnRead<-options(warn = warnDef)
    return (sa)
}



menuA<-function (choices, stInd, showLine) 
{
    ## menue function, for internal use
	nc <- length(choices)
	cat(showLine)
	if (stInd == 1) {
		op <- paste(format(seq_len(nc)), ": ", choices, sep = "")
		st <-stInd
		ed<-nc+stInd-1
	} else if (stInd == 0) {
		op <- paste(format(seq_len(nc)-1), ": ", choices, sep = "")
		st <-stInd
		ed<-nc+stInd-1
	}
	cat("  ", op, "", sep = "\n")
	repeat {
	    ind <- .Internal(menu(as.character(choices)))
	    if (ind >= st && ind <= ed) 
		return(ind)
	}
}


getinput<-function(lowlim=0,highlim=1)
{
    ## get input from user, for internal use
	repeat {
		inp<- scan(n=1, quiet = TRUE)
	    if (length(inp)) {
			if (highlim>0) {
				if (inp >= lowlim && inp <= highlim ) 
					return(inp)
			} else {
				if (inp >= lowlim )
					return(inp)
			}
	    }  
	}
}

plotRelCon<-function(BM, metaName, plotHist = FALSE, breaks, saveFig = TRUE, saveFigDir = BM$outputDir, prefixFig, rerun = FALSE)
{
    ## histogram or boxplot of relative concentration posteriors for listed metabolites with 95% credible interval
    if (missing(BM))
        return(cat("Please input batman data list.\n"))
    
    warnDef<-options("warn")$warn
    warnRead<-options(warn = -1)
    ptype = "pdf"
    n <- 4
    ns<-5
    
    m<-row.names(BM$beta)
    if (!missing(metaName))
    {
        mind<-which(!is.na(match(tolower(m),tolower(metaName)))) 
        if (length(mind)==0) 
        return(cat("No matching metabolite found...\n"))   
    } else {
        metaName <-NULL
        mind<-1:length(m)
    }
    
    ## plot batman results
    if (!is.null(BM$betaSam) && !rerun) 
    {
        r<-row.names(BM$betaSam)
        #f<-nrow(BM$betaSam)
        f<-length(mind)
        fc<-ncol(BM$betaSam)
        sno<-ncol(BM$beta)
        ind<-fc/sno
        outpdf1<-NULL
        for (i in 1:sno) 
        {
            if (plotHist)
            {
                for (j in 1:f)
                {
                    ## set subplot
                    if ((j%%n)==1)
                    {
                        if ((f-j)>=1) 
                        {
                			if (!missing(prefixFig))
                    			outpdf1 <- paste(saveFigDir, "/", prefixFig,"_spec_",i,"RelCon_", mind[j], "to",mind[min(j+n-1,length(mind))],".",ptype, sep="")
                			else
                                outpdf1 <- paste(saveFigDir,"/spec_",i,"RelCon_",mind[j], "to",mind[min(j+n-1,length(mind))],".",ptype, sep="")	
                            x11(15,7)
                			par(mfrow=c(2, 2),oma = c(0, 0, 3, 0))	
            		    } else {
                		    if (!missing(prefixFig))
                    		  	outpdf1 <- paste(saveFigDir, "/", prefixFig, "_spec_",i,"RelCon_", mind[j], ".",ptype, sep="")
                			else
                    		    outpdf1 <- paste(saveFigDir,"/spec_",i,"RelCon_",mind[j],".",ptype, sep="")
                            x11(15,7)
                            par(mfrow=c(1,1),oma = c(0, 0, 3, 0))
                        }
        	        }
        	        ## set default breaks for histogram
            		if (missing(breaks))
                    {
                        breaks <- length(BM$betaSam[mind[j],((i-1)*ind+1):(i*ind)])%/%3
                        if (breaks <=0)
                            breaks <- length(BM$betaSam[mind[j],((i-1)*ind+1):(i*ind)])
                    }    
                    hist(t(BM$betaSam[mind[j],((i-1)*ind+1):(i*ind)]),col="gray",
                    main = paste("Relative concentration\nfor ",r[mind[j]],sep=""),
                    xlab = "Relative Concentration", breaks = breaks)
            		v<-quantile(BM$betaSam[mind[j],((i-1)*ind+1):(i*ind)],p = c(2.5,97.5)/100)
                    abline(v = v, col = "red", lty=2, lwd = 2)
                    legend("topright", legend = "2.5% & 97.5% quantiles",
                    col = "red", pt.cex=2, lty =  2, cex = 0.8)
                    if (((f == j || !(j%%n))) && saveFig) 
                    {
                        title(main = list(paste("NMR Spectrum ",i, ": ", BM$specTitle[2,i], sep=""),cex=1.5,font=3), outer=TRUE)                
                        if (file.exists(outpdf1))
                            cat("Can't save figure, file", outpdf1, "already exists.\n")
                        else
                            df = dev.copy2pdf(device=x11, file = outpdf1)
                    }
                }
             } else {
                ## set subplot             
                if ((i%%n)==1)
                {
                    if ((sno-i)>=1) 
                    {
            			if (!missing(prefixFig))
                			outpdf1 <- paste(saveFigDir, "/", prefixFig,"_spec_",i,"to", min(i+n-1,sno),"RelCon.",ptype, sep="")
            			else
                            outpdf1 <- paste(saveFigDir,"/spec_",i,"to", min(i+n-1,sno),"RelCon.",ptype, sep="")	
                        x11(15,7)
            			par(mfrow=c(2, 2),oma = c(0, 0, 0, 0))	
        		    } else {
            		    if (!missing(prefixFig))
                		  	outpdf1 <- paste(saveFigDir, "/", prefixFig, "_spec_",i,"RelCon.",ptype, sep="")
            			else
                		    outpdf1 <- paste(saveFigDir,"/spec_",i,"RelCon.",ptype, sep="")
                        x11(15,7)
                        par(mfrow=c(1,1),oma = c(3, 0, 0, 0))
                    }
    	        }
                boxplot(as.data.frame(t(BM$betaSam[mind,])), 
                main = paste("NMR Spectrum ",i, ": ", BM$specTitle[2,i], sep=""), 
                names = as.character(r[mind]),
                ylab="Relative Concentration", xaxt="n") 
                axis(1,at=1:length(mind),adj=1,padj=0.5,labels=as.character(r[mind]), las=2)
             
                ## save plot
                if (((sno == i || !(i%%n))) && saveFig) 
                {
                    #title(main = list("NMR Spectrum Relative Concentration",cex=1.5,font=3), outer=TRUE)                
                    if (file.exists(outpdf1))
                        cat("Can't save figure, file", outpdf1, "already exists.\n")
                    else
                        df = dev.copy2pdf(device=x11, file = outpdf1)
                }
            }
        }
    }
    ## plot batman rerun results
    else if (!is.null(BM$betaSamRerun) && rerun) 
    {  
        #f<-nrow(BM$betaSamRerun)
        f<-length(mind)
        fc<-ncol(BM$betaSamRerun)
        sno<-length(BM$sFit)/ns
        ind<-fc/sno
        r<-row.names(BM$betaSamRerun)
        outpdf2<-NULL
        for (i in 1:sno) 
        {                                    
            if (plotHist)
            {
                for (j in 1:f)
                {   
                    ## set subplot
                    if ((j%%n)==1)
                    {
                        if ((f-j)>=1) 
                        {
                    		if (!missing(prefixFig))
                        		outpdf2 <- paste(saveFigDir, "/", prefixFig,"_RelConRerun_", mind[j], "to",mind[min(j+n-1,length(mind))],"_Rerun",".",ptype,sep="")
                    		else
                                outpdf2 <- paste(saveFigDir,"/RelConRerun_",mind[j], "to",mind[min(j+n-1,length(mind))],"_Rerun",".",ptype, sep="")	
                            x11(15,7)
                    		par(mfrow=c(2, 2),oma = c(0, 0, 3, 0))				
                        } else {
                    	    if (!missing(prefixFig))
                        		outpdf2 <- paste(saveFigDir, "/", prefixFig, "_RelConRerun_", mind[j], "_Rerun",".",ptype,sep="")
                    		else
                        	    outpdf2 <- paste(saveFigDir,"/RelConRerun_",mind[j],"_Rerun",".",ptype, sep="")
                            x11(15,7)
                    	    par(mfrow=c(1,1),oma = c(0, 0, 3, 0))	
                        }
                    }
                    ## set default breaks
            		if (missing(breaks))
                    {
                        breaks <- length(BM$betaSamRerun[mind[j],((i-1)*ind+1):(i*ind)])%/%3
                        if (breaks <=0)
                            breaks <- length(BM$betaSamRerun[mind[j],((i-1)*ind+1):(i*ind)])
                    }                
            		hist(t(BM$betaSamRerun[mind[j],((i-1)*ind+1):(i*ind)]),col="gray",
                    main = paste("Relative concentration\nfor ",r[mind[j]],sep=""),
                    xlab = "Relative Concentration", breaks = breaks)
            		v<-quantile(BM$betaSamRerun[mind[j],((i-1)*ind+1):(i*ind)],p = c(2.5,97.5)/100)
                    abline(v = v, col = "red", lty=2, lwd = 2)
                    legend("topright", legend = "2.5% & 97.5% quantiles",  col = "red",  pt.cex=2, lty =  2, cex = 0.8)
                    if (((f == j || !(j%%n))) && saveFig) 
                    {
                        title(main = list(paste("NMR Spectrum ",i, ": ", BM$specTitle[2,i], sep=""),cex=1.5,font=3), outer=TRUE)                
                        if (file.exists(outpdf2))
                            cat("Can't save figure, file", outpdf2, "already exists.\n")
                        else
                            df = dev.copy2pdf(device=x11, file = outpdf2)
                    }
                }
            } else {
                ## set subplot
                if ((i%%n)==1)
                {
                    if ((sno-i)>=1) 
                    {
            			if (!missing(prefixFig))
                			outpdf2 <- paste(saveFigDir, "/", prefixFig,"_spec_",i,"to", min(i+n-1,sno),"RelCon_Rerun.",ptype, sep="")
            			else
                            outpdf2 <- paste(saveFigDir,"/spec_",i,"to", min(i+n-1,sno),"RelCon_Rerun.",ptype,sep="")	
                        x11(15,7)
            			par(mfrow=c(2, 2),oma = c(0, 0, 3, 0))	
        		    } else {
            		    if (!missing(prefixFig))
                		  	outpdf2 <- paste(saveFigDir, "/", prefixFig, "_spec_",i,"RelCon_Rerun.",ptype, sep="")
            			else
                		    outpdf2 <- paste(saveFigDir,"/spec_",i,"RelCon_Rerun.",ptype,sep="")
                        x11(15,7)
                        par(mfrow=c(1,1),oma = c(0, 0, 3, 0))
                    }
    	        }
                boxplot(as.data.frame(t(BM$betaSamRerun[mind,])), 
                main = paste("NMR Spectrum ",i, ": ", BM$specTitle[2,i], sep=""), 
                names = as.character(r[mind]),
                ylab="Relative Concentration", xaxt="n") 
                axis(1,at=1:length(mind),adj=1,padj=0.5,labels=as.character(r[mind]), las=2)
                
                ## save plot
                if (((sno == i || !(i%%n))) && saveFig) 
                {
                    #title(main = list(paste("NMR Spectrum ",i, ": ", BM$specTitle[2,i], "(Rerun)",sep=""),cex=1.5,font=3), outer=TRUE)                
                    if (file.exists(outpdf2))
                        cat("Can't save figure, file", outpdf2, "already exists.\n")
                    else
                        df = dev.copy2pdf(device=x11, file = outpdf2)
                }
            }
        }
    } else {
        cat("No results found.\n")
   }
   warnRead<-options(warn = warnDef)
}

plotBatmanFit<-function(BM, xfrom, xto, yfrom, yto, listMeta = FALSE, metaName, saveFig = TRUE, 
                        saveFigDir = BM$outputDir, prefixFig, rerun = FALSE, placeLegend)
{
    ## plot batman metabolites fittings of NMR spectra (with down sampling)
    if (missing(BM))
    return(cat("Please input batman data list.\n"))
    
    warnDef<-options("warn")$warn
    warnRead<-options(warn = -1)
    
    if (missing(placeLegend))
    placeLegend <- "topright"
    ptype = "pdf"
    cex = 0.8
	cex1 = 0.5
    ns<-5
    ## if missing from and to parameters, use whole spectrum
    if (missing(xfrom))
    xfrom = min(BM$sFit[1,1],BM$sFit[nrow(BM$sFit),1])
    if (missing(xto))
    xto =  max(BM$sFit[1,1],BM$sFit[nrow(BM$sFit),1])
    if (xfrom > xto)
    {
        temp<-xfrom
        xfrom<-xto
        xto<-temp
    }
    if (missing(yfrom))
        yfrom <- 0
    if (missing(yto))
        ytoIP<-NULL
    pind<-which(BM$sFit[,1]<=xto & BM$sFit[,1]>=xfrom)
    if (length(pind)== 0)
        pind<-which(BM$sFit[,1]<=xfrom & BM$sFit[,1]>=xto)
    ## match metabolite 
	nometa<-FALSE
    m<-row.names(BM$beta)
    if (!missing(metaName))
    {
        mind<-which(!is.na(match(tolower(m),tolower(metaName)))) 
        if (length(mind)==0) {cat("No matching metabolite found...\n")}    
    } else if (missing(metaName) && listMeta) {
		mind<-NULL
    	nometa<-TRUE
		metaName<-"all"
 	} else {
        mind<-NULL
    	nometa<-TRUE
    	metaName<-NULL
    }

    sno<-length(BM$sFit)/ns
	n <- 2
    metaTmplty <-ns
    metaTmplwd <-2
    plotcol<- sample(rainbow(nrow(BM$beta)))
    outpdf1<-NULL
    gind<-NULL
    gapsize<-NULL
    gap<-NULL
    ## set gap for plot
    if (!is.null(BM$sFit))
    {
        x<-BM$sFit[pind,1]
        df<-abs(diff(BM$sFit[pind,1]))
        gind<-which(df>(min(df)*3))
        for (gi in 1:length(gind))
            gap<-rbind(gap, BM$sFit[pind[gind[gi]],1],BM$sFit[pind[gind[gi]+1],1])
        gap<-(sort(gap))
        lgap <- length(gap)
        gapsize<-abs(diff(gap))
        if (lgap>0)
        gapsizeodd<-gapsize[seq(1,length(gapsize),2)]
        xlim <- c(min(BM$sFit[pind,1]), max(BM$sFit[pind,1]) - gapsize[1])
        xtics <- pretty(BM$sFit[pind,1],n=(abs(diff(range(x)))%/%0.05))
        xticlab <- xtics
        lgap <- length(gap)
        lgo<-(lgap-2)/2
        if (lgap>0)  
        {
            littleones <- which(x <= gap[1])
            if (lgap > 3) 
            {
                middleones <-NULL
                for (ig in 1:lgo)
                {
                if (ig == 1)
                middleones <- list(which(x >= gap[ig*2] & x <= gap[ig*2+1]))
                else
                middleones[[ig]]<-which(x >= gap[ig*2] & x <= gap[ig*2+1])
                }                
                bigones <- which(x >= gap[lgap] )
                #lostones <- sum(c(x > gap[1] & x < gap[2],x > gap[3] & x < gap[4] ))
            } else {
                middleones <- NULL
                bigones <- which(x >= gap[2])
                #lostones <- sum(x > gap[1] & x < gap[2])
            }
            littletics <- which(xtics < gap[1])
            if (length(gapsize) > 2) 
            {
                middletics<-NULL
                maxdif<-0
                midshowat<-NULL
                midshowlab<-NULL
                middif<-NULL
                for (ig in 1:lgo)
                {
                maxdif <- maxdif + gapsizeodd[ig]
                middif<-rbind(middif, maxdif)
                middletics <- which(xtics >= gap[ig*2] & xtics <=gap[ig*2+1])
                midshowat<-c(midshowat, xtics[middletics] - maxdif)
                midshowlab<-c(midshowlab,xticlab[middletics])
                }
                maxdif<-maxdif + gapsizeodd[length(gapsizeodd)]
                xlim <- c(min(x), max(x)-maxdif)
                bigtics <- which(xtics >= gap[lgap])
                show.at <- c(xtics[littletics], midshowat, xtics[bigtics] - maxdif)
                #show.at <- c(xtics[littletics], xtics[middletics] - gapsize[1], xtics[bigtics] - (gapsize[1] + gapsize[3]))
                show.labels <- c(xticlab[littletics], midshowlab,xticlab[bigtics])
            } else {
                xlim <- c(min(x), max(x) - gapsize[1])
                bigtics <- which(xtics >= gap[2])
                show.at <- c(xtics[littletics], xtics[bigtics] - gapsize[1])
                show.labels <- c(xticlab[littletics], xticlab[bigtics])
            }
        }
    }
    ## plot batman results
    if (!is.null(BM$sFit) && !rerun) 
    {
    	for (j in 1:sno)
        {
    		if ((j%%n) == 1)
            { 
                if ((sno-j)>=1) 
                {
        			if (!missing(prefixFig))
            			outpdf1 <- paste(saveFigDir, "/", prefixFig,"_specFit_", j, "to",j+n-1,"_",metaName,".",ptype, sep="")
        			else
                        outpdf1 <- paste(saveFigDir,"/specFit_",j, "to",j+n-1,"_",metaName,".",ptype, sep="")	           
                    x11(15,7)
                    par(mfrow=c(n,1))	         		
    		    } else {
        		    if (!missing(prefixFig))
            			outpdf1 <- paste(saveFigDir, "/", prefixFig, "_specFit_", j, "_",metaName,".",ptype, sep="")
        			else
            		    outpdf1 <- paste(saveFigDir,"/specFit_",j,"_",metaName,".",ptype, sep="")
                    x11(15,7)
                }
            } 
    		i <- (ns*(j-1)+1)
    
            if (is.null(ytoIP))
                yto <- max(max(BM$sFit[pind,i+1]),max(BM$sFit[pind,i+2]),max(BM$sFit[pind,i+3]),max(BM$sFit[pind,i+4]))
    
    		if (yfrom > yto)
            {
                temp<-yfrom
                yfrom<-yto
                yto<-temp
            }      
            if (length(gap)>0)
            {    
                ## plot metabolites fit with gap
                ylim <- c(yfrom, yto)
                ytics <- pretty(ylim)
                yticlab <- ytics
                
                plot(BM$sFit[pind[littleones],i], BM$sFit[pind[littleones],i+1], xlim = rev(xlim), 
                ylim = ylim, axes = FALSE, lwd = 0.5, col = 4, lty = 1,
                type="l", xlab="ppm", ylab="Standardized Intensity",
                main=paste("NMR Spectrum ",j,": ",BM$specTitle[2,j],sep=""))
                box()
                axis(2, at = ytics, labels = yticlab)
                axis(1, at = show.at, labels = show.labels)
                axis.break(1, gap[1], style = "zigzag")
                if (length(gapsize) > 2) {
                    for (ig in 1:lgo)
                    {
                    axis.break(1, gap[ig*2+1] - middif[ig], style = "zigzag")
                    lines(BM$sFit[pind[middleones[[ig]]],i]-middif[ig],BM$sFit[pind[middleones[[ig]]],i+1],col=4,lwd=0.5,lty=1)
                    lines(BM$sFit[pind[middleones[[ig]]],i]-middif[ig],BM$sFit[pind[middleones[[ig]]],i+2],col=3,lwd=0.5,lty=1)
                    lines(BM$sFit[pind[middleones[[ig]]],i]-middif[ig],BM$sFit[pind[middleones[[ig]]],i+3],col=2,lwd=0.5,lty=1)
                    lines(BM$sFit[pind[middleones[[ig]]],i]-middif[ig],BM$sFit[pind[middleones[[ig]]],i+4],col=1,lwd=0.5,lty=1)
                    }
                    lines(BM$sFit[pind[bigones],i]-maxdif,BM$sFit[pind[bigones],i+1],col=4,lwd=0.5,lty=1)
                    lines(BM$sFit[pind[littleones],i],BM$sFit[pind[littleones],i+2],col=3,lwd=0.5,lty=1)
                    lines(BM$sFit[pind[bigones],i]-maxdif,BM$sFit[pind[bigones],i+2],col=3,lwd=0.5,lty=1)
                    lines(BM$sFit[pind[littleones],i],BM$sFit[pind[littleones],i+3],col=2,lwd=0.5,lty=1)
                    lines(BM$sFit[pind[bigones],i]-maxdif,BM$sFit[pind[bigones],i+3],col=2,lwd=0.5,lty=1)
                    lines(BM$sFit[pind[littleones],i],BM$sFit[pind[littleones],i+4],col=1,lwd=0.5,lty=1)
                    lines(BM$sFit[pind[bigones],i]-maxdif,BM$sFit[pind[bigones],i+4],col=1,lwd=0.5,lty=1)
                    if (listMeta && nometa)
            	    {
                        ## if listMeta is TRUE and missing metabolite name, plot all
                        for (i2 in 1:nrow(BM$beta))
                        {
                            ytmp <- BM$beta[i2,j]*BM$metaTemp[pind,i2+(j-1)*nrow(BM$beta)]
                            lines(BM$sFit[pind[littleones],i],ytmp[littleones],col=plotcol[i2],lwd = metaTmplwd, lty = metaTmplty)
                            for (ig in 1:lgo)
                            lines(BM$sFit[pind[middleones[[ig]]],i]-middif[ig],ytmp[middleones[[ig]]],col=plotcol[i2],lwd = metaTmplwd, lty = metaTmplty)
                            
                            lines(BM$sFit[pind[bigones],i]-maxdif,ytmp[bigones],col=plotcol[i2],lwd = metaTmplwd, lty = metaTmplty)
                        }
                        legend(placeLegend, c("Original Spectrum", "Metabolites Fit", "Wavelet Fit", "Fit Sum",
                        row.names(BM$beta)), col=c(4,3,2,1,plotcol), ncol = 2, cex = cex1,
                        lty=c(1,1,1,1,rep(metaTmplty, nrow(BM$beta))),lwd = c(0.5,0.5,0.5,0.5,rep(metaTmplwd, nrow(BM$beta))))
                    } else if (length(mind)!=0) {
                        ## plot named metabolite
                        ## lines(BM$beta[mind,j]*BM$metaTemp[pind,mind+(j-1)*nrow(BM$beta)],col=plotcol[mind], lwd = metaTmplwd, lty = metaTmplty )
                        ytmp <- BM$beta[mind,j]*BM$metaTemp[pind,mind+(j-1)*nrow(BM$beta)]
                        lines(BM$sFit[pind[littleones],i],ytmp[littleones],col=plotcol[mind], lwd = metaTmplwd, lty = metaTmplty)
                        for (ig in 1:lgo)
                        lines(BM$sFit[pind[middleones[[ig]]],i]-middif[ig],ytmp[middleones[[ig]]],col=plotcol[mind], lwd = metaTmplwd, lty = metaTmplty)
                        lines(BM$sFit[pind[bigones],i]-maxdif,ytmp[bigones],col=plotcol[mind], lwd = metaTmplwd, lty = metaTmplty)
                                        
                        legend(placeLegend, c("Original Spectrum", "Metabolites Fit", "Wavelet Fit", "Fit Sum",
                        row.names(BM$beta)[mind]), col=c(4,3,2,1,plotcol[mind]), cex = cex,
                        lty=c(1,1,1,1,rep(metaTmplty, length(mind))),lwd = c(0.5,0.5,0.5,0.5,rep(metaTmplwd, length(mind))))
                    } else {
                        legend(placeLegend, c("Original Spectrum", "Metabolites Fit", "Wavelet Fit", "Fit Sum"), col=c(4,3,2,1), 
                        lty=c(1,1,1,1),lwd = c(0.5,0.5,0.5,0.5), cex = cex)
                    }
                } else {
                    lines(BM$sFit[pind[bigones],i]-(gapsize[1]),BM$sFit[pind[bigones],i+1],col=4,lwd=0.5,lty=1)
                    lines(BM$sFit[pind[littleones],i],BM$sFit[pind[littleones],i+2],col=3,lwd=0.5,lty=1)
                    lines(BM$sFit[pind[bigones],i]-(gapsize[1]),BM$sFit[pind[bigones],i+2],col=3,lwd=0.5,lty=1)
                    lines(BM$sFit[pind[littleones],i],BM$sFit[pind[littleones],i+3],col=2,lwd=0.5,lty=1)
                    lines(BM$sFit[pind[bigones],i]-(gapsize[1]),BM$sFit[pind[bigones],i+3],col=2,lwd=0.5,lty=1)
                    lines(BM$sFit[pind[littleones],i],BM$sFit[pind[littleones],i+4],col=1,lwd=0.5,lty=1)
                    lines(BM$sFit[pind[bigones],i]-(gapsize[1]),BM$sFit[pind[bigones],i+4],col=1,lwd=0.5,lty=1)

                    if (listMeta && nometa)
            	    {
                        ## if listMeta is TRUE and missing metabolite name, plot all
                        for (i2 in 1:nrow(BM$beta))
                        {
                            ytmp <- BM$beta[i2,j]*BM$metaTemp[pind,i2+(j-1)*nrow(BM$beta)]
                            lines(BM$sFit[pind[littleones],i],ytmp[littleones],col=plotcol[i2],lwd = metaTmplwd, lty = metaTmplty)
                            lines(BM$sFit[pind[bigones],i]-(gapsize[1]),ytmp[bigones],col=plotcol[i2],lwd = metaTmplwd, lty = metaTmplty)
                        }
                        legend(placeLegend, c("Original Spectrum", "Metabolites Fit", "Wavelet Fit", "Fit Sum",
                        row.names(BM$beta)), col=c(4,3,2,1,plotcol), ncol = 2, cex = cex1,
                        lty=c(1,1,1,1,rep(metaTmplty, nrow(BM$beta))),lwd = c(0.5,0.5,0.5,0.5,rep(metaTmplwd, nrow(BM$beta))))
                    } else if (length(mind)!=0) {
                        ## plot named metabolite
                        ## lines(BM$beta[mind,j]*BM$metaTemp[pind,mind+(j-1)*nrow(BM$beta)],col=plotcol[mind], lwd = metaTmplwd, lty = metaTmplty )
                        ytmp <- BM$beta[mind,j]*BM$metaTemp[pind,mind+(j-1)*nrow(BM$beta)]
                        lines(BM$sFit[pind[littleones],i],ytmp[littleones],col=plotcol[mind], lwd = metaTmplwd, lty = metaTmplty)
                        lines(BM$sFit[pind[bigones],i]-(gapsize[1]),ytmp[bigones],col=plotcol[mind], lwd = metaTmplwd, lty = metaTmplty)
                                        
                        legend(placeLegend, c("Original Spectrum", "Metabolites Fit", "Wavelet Fit", "Fit Sum",
                        row.names(BM$beta)[mind]), col=c(4,3,2,1,plotcol[mind]), cex = cex,
                        lty=c(1,1,1,1,rep(metaTmplty, length(mind))),lwd = c(0.5,0.5,0.5,0.5,rep(metaTmplwd, length(mind))))
                    } else {
                        legend(placeLegend, c("Original Spectrum", "Metabolites Fit", "Wavelet Fit", "Fit Sum"), col=c(4,3,2,1), 
                            lty=c(1,1,1,1),lwd = c(0.5,0.5,0.5,0.5), cex = cex)
                    }
                }
            } else {
                ## plot metabolites fit without gap
                plot(BM$sFit[pind,i],BM$sFit[pind,i+1],type="l",xlim=rev(range(BM$sFit[pind,i])),xlab="ppm",
                ylab="Standardized Intensity", main=paste("NMR Spectrum ",j,": ",BM$specTitle[2,j],sep=""), 
                ylim = c(yfrom, yto), lwd = 0.5, col = 4, lty = 1)
                lines(BM$sFit[pind,i],BM$sFit[pind,i+2],col=3, lwd = 0.5, lty = 1)
        		lines(BM$sFit[pind,i],BM$sFit[pind,i+3],col=2, lwd = 0.5, lty = 1)
        		lines(BM$sFit[pind,i],BM$sFit[pind,i+4],col=1, lwd = 0.5, lty = 1)
    
            	if (listMeta && nometa)
            	{
                    ## if listMeta is TRUE and missing metabolite name, plot all
                    for (i2 in 1:nrow(BM$beta))
                    lines(BM$sFit[pind,i],BM$beta[i2,j]*BM$metaTemp[pind,i2+(j-1)*nrow(BM$beta)],col=plotcol[i2], lwd = metaTmplwd, lty = metaTmplty )
                    
                    legend(placeLegend, c("Original Spectrum", "Metabolites Fit", "Wavelet Fit", "Fit Sum",
                    row.names(BM$beta)), col=c(4,3,2,1,plotcol), ncol = 2, cex = cex1,
                    lty=c(1,1,1,1,rep(metaTmplty, nrow(BM$beta))),lwd = c(0.5,0.5,0.5,0.5,rep(metaTmplwd, nrow(BM$beta))))
                } else if (length(mind)!=0) {
                    ## plot named metabolite
                    lines(BM$sFit[pind,i],BM$beta[mind,j]*BM$metaTemp[pind,mind+(j-1)*nrow(BM$beta)],col=plotcol[mind], lwd = metaTmplwd, lty = metaTmplty )
                    
                    legend(placeLegend, c("Original Spectrum", "Metabolites Fit", "Wavelet Fit", "Fit Sum",
                    row.names(BM$beta)[mind]), col=c(4,3,2,1,plotcol[mind]), cex = cex,
                    lty=c(1,1,1,1,rep(metaTmplty, length(mind))),lwd = c(0.5,0.5,0.5,0.5,rep(metaTmplwd, length(mind))))
                } else {
                    legend(placeLegend, c("Original Spectrum", "Metabolites Fit", "Wavelet Fit", "Fit Sum"), col=c(4,3,2,1), 
                    lty=c(1,1,1,1),lwd = c(0.5,0.5,0.5,0.5), cex = cex)
                }
            }
            if ((sno == j || !(j%%n)) && saveFig) {
                if (file.exists(outpdf1))
                    cat("Can't save figure, file", outpdf1, "already exists.\n")
                else
                    df = dev.copy2pdf(device=x11, file = outpdf1)
            }
        }
    }
    ## plot batman rerun results
	else if (!is.null(BM$sFitRerun) && rerun) 
    {     
        outpdf2 <-NULL
		for (j in 1:sno)
        {
            ## set subplot
			if ((j%%n)==1)
            {
                if ((sno-j)>=1) 
                {
        			if (!missing(prefixFig))
            			outpdf2 <- paste(saveFigDir, "/", prefixFig,"_specfitRerun_", j, "to",j+n-1,"_",metaName,".",ptype, sep="")
        			else
                        outpdf2 <- paste(saveFigDir,"/specfitRerun_",j, "to",j+n-1,"_",metaName,".",ptype, sep="")
                    x11(15,7)
        			par(mfrow=c(n,1))			
    			} else {
                    x11(15,7)
        			if (!missing(prefixFig))
            			outpdf2 <- paste(saveFigDir, "/", prefixFig,"_specfitRerun_", j, "_",metaName,".",ptype, sep="")
        			else
            		    outpdf2 <- paste(saveFigDir,"/specfitRerun_",j,"_",metaName,".",ptype, sep="")
                } 
            }
			i = ns*(j-1)+1

            if (is.null(ytoIP))
                yto <- max(max(BM$sFitRerun[pind,i+1]),max(BM$sFitRerun[pind,i+2]),max(BM$sFitRerun[pind,i+3]),max(BM$sFitRerun[pind,i+4]))
            
    		if (yfrom > yto)
            {
                temp<-yfrom
                yfrom<-yto
                yto<-temp
            }
            if (length(gap)>0)
            {    
                ## plot metabolites fit with gap
                ylim <- c(yfrom, yto)
                ytics <- pretty(ylim)
                yticlab <- ytics
                
                plot(BM$sFitRerun[pind[littleones],i], BM$sFitRerun[pind[littleones],i+1], xlim = rev(xlim), 
                ylim = ylim, axes = FALSE, lwd = 0.5, col = 4, lty = 1,type="l", xlab="ppm", ylab="Standardized Intensity",
                main=paste("NMR Spectrum ",j, ": ",BM$specTitle[2,j],"(Rerun)", sep=""))
                box()
                axis(2, at = ytics, labels = yticlab)
                axis(1, at = show.at, labels = show.labels)
                axis.break(1, gap[1], style = "zigzag")
                if (length(gapsize) > 2) {
                    for (ig in 1:lgo)
			        {
                    axis.break(1, gap[ig*2+1] - middif[ig], style = "zigzag")
                    lines(BM$sFitRerun[pind[middleones[[ig]]],i]-middif[ig],BM$sFitRerun[pind[middleones[[ig]]],i+1],col=4,lwd=0.5,lty=1)
                    lines(BM$sFitRerun[pind[middleones[[ig]]],i]-middif[ig],BM$sFitRerun[pind[middleones[[ig]]],i+2],col=3,lwd=0.5,lty=1)
                    lines(BM$sFitRerun[pind[middleones[[ig]]],i]-middif[ig],BM$sFitRerun[pind[middleones[[ig]]],i+3],col=2,lwd=0.5,lty=1)
                    lines(BM$sFitRerun[pind[middleones[[ig]]],i]-middif[ig],BM$sFitRerun[pind[middleones[[ig]]],i+4],col=1,lwd=0.5,lty=1)
                    }
                    lines(BM$sFitRerun[pind[bigones],i]-maxdif,BM$sFitRerun[pind[bigones],i+1],col=4,lwd=0.5,lty=1)
                    lines(BM$sFitRerun[pind[littleones],i],BM$sFitRerun[pind[littleones],i+2],col=3,lwd=0.5,lty=1)
                    lines(BM$sFitRerun[pind[bigones],i]-maxdif,BM$sFitRerun[pind[bigones],i+2],col=3,lwd=0.5,lty=1)
                    lines(BM$sFitRerun[pind[littleones],i],BM$sFitRerun[pind[littleones],i+3],col=2,lwd=0.5,lty=1)
                    lines(BM$sFitRerun[pind[bigones],i]-maxdif,BM$sFitRerun[pind[bigones],i+3],col=2,lwd=0.5,lty=1)
                    lines(BM$sFitRerun[pind[littleones],i],BM$sFitRerun[pind[littleones],i+4],col=1,lwd=0.5,lty=1)
                    lines(BM$sFitRerun[pind[bigones],i]-maxdif,BM$sFitRerun[pind[bigones],i+4],col=1,lwd=0.5,lty=1)
                    
                    if (listMeta && nometa)
            	    {
                        ## if listMeta is TRUE and missing metabolite name, plot all
                        for (i2 in 1:nrow(BM$betaRerun))
                        {
                            ytmp <- BM$betaRerun[i2,j]*BM$metaTempRerun[pind,i2+(j-1)*nrow(BM$betaRerun)]
                            lines(BM$sFitRerun[pind[littleones],i],ytmp[littleones],col=plotcol[i2],lwd = metaTmplwd, lty = metaTmplty)
                            for (ig in 1:lgo)
                            lines(BM$sFitRerun[pind[middleones[[ig]]],i]-middif[ig],ytmp[middleones[[ig]]],col=plotcol[i2],lwd = metaTmplwd, lty = metaTmplty)
                            lines(BM$sFitRerun[pind[bigones],i]-maxdif,ytmp[bigones],col=plotcol[i2],lwd = metaTmplwd, lty = metaTmplty)
                        }
                        legend(placeLegend, c("Original Spectrum", "Metabolites Fit", "Wavelet Fit", "Fit Sum",
                        row.names(BM$betaRerun)), col=c(4,3,2,1,plotcol), ncol = 2, cex = cex1,
                        lty=c(1,1,1,1,rep(metaTmplty, nrow(BM$betaRerun))),lwd = c(0.5,0.5,0.5,0.5,rep(metaTmplwd, nrow(BM$betaRerun))))
                    } else if (length(mind)!=0) {
                        ## plot named metabolite
                        ## lines(BM$beta[mind,j]*BM$metaTemp[pind,mind+(j-1)*nrow(BM$beta)],col=plotcol[mind], lwd = metaTmplwd, lty = metaTmplty )
                        ytmp <- BM$betaRerun[mind,j]*BM$metaTempRerun[pind,mind+(j-1)*nrow(BM$betaRerun)]
                        lines(BM$sFitRerun[pind[littleones],i],ytmp[littleones],col=plotcol[mind], lwd = metaTmplwd, lty = metaTmplty)
                        for (ig in 1:lgo)
                        lines(BM$sFitRerun[pind[middleones[[ig]]],i]-middif[ig],ytmp[middleones[[ig]]],col=plotcol[mind], lwd = metaTmplwd, lty = metaTmplty)
                        lines(BM$sFitRerun[pind[bigones],i]-maxdif,ytmp[bigones],col=plotcol[mind], lwd = metaTmplwd, lty = metaTmplty)
                                        
                        legend(placeLegend, c("Original Spectrum", "Metabolites Fit", "Wavelet Fit", "Fit Sum",
                        row.names(BM$betaRerun)[mind]), col=c(4,3,2,1,plotcol[mind]), cex = cex,
                        lty=c(1,1,1,1,rep(metaTmplty, length(mind))),lwd = c(0.5,0.5,0.5,0.5,rep(metaTmplwd, length(mind))))
                    } else {
                        legend(placeLegend, c("Original Spectrum", "Metabolites Fit", "Wavelet Fit", "Fit Sum"), col=c(4,3,2,1), 
                        lty=c(1,1,1,1),lwd = c(0.5,0.5,0.5,0.5), cex = cex)
                    }
                } else {
                    lines(BM$sFitRerun[pind[bigones],i]-(gapsize[1]),BM$sFitRerun[pind[bigones],i+1],col=4,lwd=0.5,lty=1)
                    lines(BM$sFitRerun[pind[littleones],i],BM$sFitRerun[pind[littleones],i+2],col=3,lwd=0.5,lty=1)
                    lines(BM$sFitRerun[pind[bigones],i]-(gapsize[1]),BM$sFitRerun[pind[bigones],i+2],col=3,lwd=0.5,lty=1)
                    lines(BM$sFitRerun[pind[littleones],i],BM$sFitRerun[pind[littleones],i+3],col=2,lwd=0.5,lty=1)
                    lines(BM$sFitRerun[pind[bigones],i]-(gapsize[1]),BM$sFitRerun[pind[bigones],i+3],col=2,lwd=0.5,lty=1)
                    lines(BM$sFitRerun[pind[littleones],i],BM$sFitRerun[pind[littleones],i+4],col=1,lwd=0.5,lty=1)
                    lines(BM$sFitRerun[pind[bigones],i]-(gapsize[1]),BM$sFitRerun[pind[bigones],i+4],col=1,lwd=0.5,lty=1)

                    if (listMeta && nometa)
            	    {
                        ## if listMeta is TRUE and missing metabolite name, plot all
                        for (i2 in 1:nrow(BM$betaRerun))
                        {
                            ytmp <- BM$betaRerun[i2,j]*BM$metaTempRerun[pind,i2+(j-1)*nrow(BM$betaRerun)]
                            lines(BM$sFitRerun[pind[littleones],i],ytmp[littleones],col=plotcol[i2],lwd = metaTmplwd, lty = metaTmplty)
                            lines(BM$sFitRerun[pind[bigones],i]-(gapsize[1]),ytmp[bigones],col=plotcol[i2],lwd = metaTmplwd, lty = metaTmplty)
                        }
                        legend(placeLegend, c("Original Spectrum", "Metabolites Fit", "Wavelet Fit", "Fit Sum",
                        row.names(BM$betaRerun)), col=c(4,3,2,1,plotcol), ncol = 2, cex = cex1,
                        lty=c(1,1,1,1,rep(metaTmplty, nrow(BM$betaRerun))),lwd = c(0.5,0.5,0.5,0.5,rep(metaTmplwd, nrow(BM$betaRerun))))
                    } else if (length(mind)!=0) {
                        ## plot named metabolite
                        ## lines(BM$beta[mind,j]*BM$metaTemp[pind,mind+(j-1)*nrow(BM$beta)],col=plotcol[mind], lwd = metaTmplwd, lty = metaTmplty )
                        ytmp <- BM$betaRerun[mind,j]*BM$metaTempRerun[pind,mind+(j-1)*nrow(BM$betaRerun)]
                        lines(BM$sFitRerun[pind[littleones],i],ytmp[littleones],col=plotcol[mind], lwd = metaTmplwd, lty = metaTmplty)
                        lines(BM$sFitRerun[pind[bigones],i]-(gapsize[1]),ytmp[bigones],col=plotcol[mind], lwd = metaTmplwd, lty = metaTmplty)
                                        
                        legend(placeLegend, c("Original Spectrum", "Metabolites Fit", "Wavelet Fit", "Fit Sum",
                        row.names(BM$betaRerun)[mind]), col=c(4,3,2,1,plotcol[mind]), cex = cex,
                        lty=c(1,1,1,1,rep(metaTmplty, length(mind))),lwd = c(0.5,0.5,0.5,0.5,rep(metaTmplwd, length(mind))))
                    } else {
                        legend(placeLegend, c("Original Spectrum", "Metabolites Fit", "Wavelet Fit", "Fit Sum"), col=c(4,3,2,1), 
                        lty=c(1,1,1,1),lwd = c(0.5,0.5,0.5,0.5), cex = cex)
                    }
                }
            } else {
				## plot metabolites fit without gap
        		plot(BM$sFitRerun[pind,i],BM$sFitRerun[pind,i+1],type="l",xlim=rev(range(BM$sFitRerun[pind,i])),xlab="ppm",
        			 ylab="Standardized Intensity", main=paste("NMR Spectrum ",j, ": ",BM$specTitle[2,j],"(Rerun)", sep=""), 
                     ylim = c(yfrom, yto), lwd = 0.5, col = 4, lty = 1)
                #axis(1, 1:length(pind), lab = format(BM$sFitRerun[pind,i]),xlim=rev(range(BM$sFitRerun[pind,i])))     
        		lines(BM$sFitRerun[pind,i],BM$sFitRerun[pind,i+2],col=3, lwd = 0.5, lty = 1)
        		lines(BM$sFitRerun[pind,i],BM$sFitRerun[pind,i+3],col=2, lwd = 0.5, lty = 1)
        		lines(BM$sFitRerun[pind,i],BM$sFitRerun[pind,i+4],col=1, lwd = 0.5, lty = 1)
                
        	    if (listMeta  && nometa)
        	    {
                    ## if listMeta is TRUE and missing metabolite name, plot all
                    for (i2 in 1:nrow(BM$betaRerun))
                    lines(BM$sFitRerun[pind,i],BM$betaRerun[i2,j]*BM$metaTempRerun[pind,i2+(j-1)*nrow(BM$betaRerun)],col=plotcol[i2], lwd = metaTmplwd, lty = metaTmplty )
                    legend(placeLegend, c("Original Spectrum", "Metabolites Fit", "Wavelet Fit", "Fit Sum",
                    row.names(BM$betaRerun)), col=c(4,3,2,1,plotcol), ncol = 2, cex = cex1,
                    lty=c(1,1,1,1,rep(metaTmplty, nrow(BM$betaRerun))),lwd = c(0.5,0.5,0.5, 0.5,rep(metaTmplwd, nrow(BM$betaRerun))))
                } else if (length(mind)!=0) {
                    ## plot named metabolite
                    lines(BM$sFitRerun[pind,i],BM$betaRerun[mind,j]*BM$metaTempRerun[pind,mind+(j-1)*nrow(BM$betaRerun)],col=plotcol[mind], lwd = metaTmplwd, lty = metaTmplty )
                    legend(placeLegend, c("Original Spectrum", "Metabolites Fit", "Wavelet Fit", "Fit Sum",
                    row.names(BM$betaRerun)[mind]), col=c(4,3,2,1,plotcol[mind]), cex = cex,
                    lty=c(1,1,1,1,rep(metaTmplty, length(mind))),lwd = c(0.5,0.5,0.5,0.5,rep(metaTmplwd, length(mind))))
                } else {
                    legend(placeLegend, c("Original Spectrum", "Metabolites Fit", "Wavelet Fit", "Fit Sum"), col=c(4,3,2,1), 
                    lty=c(1,1,1,1),lwd = c(0.5,0.5,0.5,0.5), cex = cex)
                }
            }
            ## save plot
            if ((sno == j || !(j%%n)) && saveFig) 
            {
                if (file.exists(outpdf2))
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

newDir<-function(runBATMANDir = getwd(), overwriteFile = FALSE)
{
    ## create a new directory with installed data files
	dirRoot<-paste(runBATMANDir,"/runBATMAN",sep="")
	if(!file.exists(dirRoot)) {
		dir.create(dirRoot)
	}
	dirIP<-paste(dirRoot,"/BatmanInput",sep="")
	if(!file.exists(dirIP)){
		dir.create(dirIP)
    }
	dirOP<-paste(dirRoot,"/BatmanOutput",sep="")
	if(!file.exists(dirOP)) {
		dir.create(dirOP)
	}
	extdir<-system.file("extdata",package="batman")
	dirbo<-paste(extdir,"/batmanOptions.txt",sep="")
	cp<-file.copy(dirbo,dirIP, overwrite=overwriteFile)
	dird<-paste(extdir,"/NMRdata.txt",sep="")
	cp<-file.copy(dird,dirIP, overwrite=overwriteFile)
	dirL<-paste(extdir,"/metabolitesList.csv",sep="")
	cp<-file.copy(dirL,dirIP, overwrite=overwriteFile)
	dirR<-paste(extdir,"/multi_data.csv",sep="")
	cp<-file.copy(dirR,dirIP,overwrite=overwriteFile)
	dirRU<-paste(extdir,"/multi_data_user.csv",sep="")
	cp<-file.copy(dirRU,dirIP,overwrite=overwriteFile)
	
	dirA<-c(dirRoot, dirIP, dirOP)
	return (dirA)
}

readBatmanOutput<-function(dirOP) 
{
    warnDef<-options("warn")$warn
    warnRead<-options(warn = -1)
    ## reads in batman output data files 
    con  <- file(paste(dirOP,"/batmanOptions.txt",sep=""), open = "r")
    oneLine <- readLines(con, n = 30, warn = FALSE)
    fL<-substr(oneLine,1,1)
    nL<-which(is.na(match(fL,"%")))
    myVector <- strsplit(oneLine[nL[2]], ":")
    NoSpectra <- as.numeric(myVector[[1]][2])
    close(con)
    
    specTitle<-read.table(paste(dirOP,"/spectraTitle.txt",sep=""), header=FALSE,sep="\t")
    mL<-read.csv(paste(dirOP,"/metabolitesListUsed.txt",sep=""), header=F,colClasses="character")
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
 
    fdir <- paste(dirOP, "/specFit_1_rr_0.txt", sep="")
    if (file.exists(fdir))
    {                
        if (!file.info(fdir)$size == 0) 
        {
            r<-read.table(fdir,sep = "\t",header=T)
        } 
    }
    fdir <- paste(dirOP, "/NMRdata_mod_1.txt", sep="")
    if (file.exists(fdir))
    {                
        if (!file.info(fdir)$size == 0) 
        {
            Ndata<-read.table(fdir, header = FALSE, sep = "\t")
        }   
    } 
    fdir <- paste(dirOP, "/specFitHR_1_rr_0.txt", sep="")
    if (file.exists(fdir))
    {                
        if (!file.info(fdir)$size == 0) 
        {
            rH<-read.table(fdir,sep = "\t",header=F)
            rH<-cbind(Ndata,rH)
        }
    }
    fdir <- paste(dirOP, "/metaTempHR_1_rr_0.txt", sep="")
    if (file.exists(fdir))
    {                
        if (!file.info(fdir)$size == 0) 
        {
            LH<-read.table(fdir,sep = "\t", header=F)
            if (length(mL) == ncol(LH))
                names(LH)<-mL
            if (nrow(mL) == ncol(LH))
                names(LH)<-t(mL)
        }
    }
    fdir <- paste(dirOP, "/metaTemp_1_rr_0.txt", sep="")
    if (file.exists(fdir))
    {                
        if (!file.info(fdir)$size == 0) 
        {
            L<-read.table(fdir,sep = "\t", header=F)
            if (length(mL) == ncol(L))
                names(L)<-mL
            if (nrow(mL) == ncol(L))
                names(L)<-t(mL)
        }
    }
    fdir <- paste(dirOP, "/beta_1_rr_0.txt", sep="")
    if (file.exists(fdir))
    {                
        if (!file.info(fdir)$size == 0) 
        {
            bet<-read.table(fdir,header=F)
            if (length(mL) == nrow(bet))
                row.names(bet)<-mL
            if ( nrow(mL) == nrow(bet) )
                row.names(bet)<-t(mL)
        }
    }
    fdir <- paste(dirOP, "/theta_sam_1_rr_0.txt", sep="") 
    if (file.exists(fdir))
    {                
        if (!file.info(fdir)$size == 0) 
        {
            if (!(file.info(fdir)$size == 0))
            {
            	the<-read.table(fdir,sep = "\t",header=F) 
            }
        }
    }
    fdir <- paste(dirOP, "/metaFit_sam_1_rr_0.txt", sep="") 
    if (file.exists(fdir))
    {                
        if (!file.info(fdir)$size == 0) 
        {
        	sfitsam<-read.table(fdir,sep = "\t",header=F)
        }  
    } 
    fdir <- paste(dirOP, "/lambda_sam_1_rr_0.txt", sep="")
    if (file.exists(fdir))
    {                
        if (!file.info(fdir)$size == 0) 
        {
            lam<-read.table(fdir,header=F)
        }    
    }   
    fdir <- paste(dirOP, "/beta_sam_1_rr_0.txt", sep="")
    if (file.exists(fdir))
    {                
        if (!file.info(fdir)$size == 0) 
        {
            betasam<-read.table(fdir,sep = "\t",header=F)
            if (length(mL) == nrow(betasam))
                row.names(betasam)<-mL
            if (nrow(mL) == nrow(betasam))
                row.names(betasam)<-t(mL)
        }  
    }
    fdir <- paste(dirOP, "/delta_sam_1_rr_0.txt", sep="")
    if (file.exists(fdir))
    {             
        if (!file.info(fdir)$size == 0) 
        {
            del<-read.delim(fdir,header=T)
        }
    }
    fdir <- paste(dirOP, "/delta_sam_1.txt", sep="")
    if (file.exists(fdir))
    {             
        if (!file.info(fdir)$size == 0) 
        {
            del<-read.delim(fdir,header=T)
        }
    }
    fdir <- paste(dirOP, "/delta_draw_mean_1.txt", sep="")
    if (file.exists(fdir))
    {                
        if (!file.info(fdir)$size == 0) 
        {
            delmean<-read.delim(fdir,header=F)
            if (ncol(del) == nrow(delmean))
                row.names(delmean)<-names(del)
        }
    }
    ## read in individual metabolite fit posteriors
    brow<-nrow(betasam)
    bcol<-ncol(betasam)
    fdir <- paste(dirOP, "/metaIndFit_sam_1_rr_0.txt", sep="")
    if (file.exists(fdir))
    { 
        if (!file.info(fdir)$size == 0)
        {
             metafitsamInd<-read.table(fdir,sep="\t",header=F)  
             for (j in 1:bcol)
             {
                 if (j == 1)
                     metaindfitsam<-t(t(metafitsamInd[,((j-1)*brow+1):(j*brow)])*betasam[,j])
                 else
                     metaindfitsam<-cbind(metaindfitsam, t(t(metafitsamInd[,((j-1)*brow+1):(j*brow)])*betasam[,j]))
             }                  
        }
    }
    ## more than 1 spectra
    if (NoSpectra>1)
    {
        for (i in 2:NoSpectra)   
        { 
            fdir <- paste(dirOP, "/specFit_", i,"_rr_0.txt", sep="")
            if (file.exists(fdir))
            {
                if (!file.info(fdir)$size == 0)
                {
                    r<-cbind(r, read.table(fdir,sep = "\t",header=T))
                }
            }    
            fdir <- paste(dirOP, "/NMRdata_mod_",i,".txt", sep="")
            if (file.exists(fdir))
            {
                if (!file.info(fdir)$size == 0)
                {
                    Ndata<-read.table(fdir, header = FALSE, sep = "\t")
                }
            }
            fdir <- paste(dirOP, "/specFitHR_", i,"_rr_0.txt", sep="")
            if (file.exists(fdir))
            { 
                if (!file.info(fdir)$size == 0)
                {               
                    rH<-cbind(rH, Ndata, read.table(fdir,sep = "\t",header=F))
                }
            }
            fdir <- paste(dirOP,  "/metaTempHR_",i, "_rr_0.txt", sep="")
            if (file.exists(fdir))
            {
                if (!file.info(fdir)$size == 0)
                {
                    LtH<-read.table(fdir,sep = "\t",header=F)
                    if (length(mL) == ncol(LtH))
                        names(LtH)<-mL 
                    if (nrow(mL) == ncol(LtH))
                        names(LtH)<-t(mL) 
                    LH<-cbind(LH, LtH)
                }
            }            
            fdir <- paste(dirOP,  "/beta_",i, "_rr_0.txt", sep="")
            if (file.exists(fdir))
            {
                if (!file.info(fdir)$size == 0)
                {
                    bet<-cbind(bet, read.table(fdir,header=F))
                }
            }
            fdir <- paste(dirOP,  "/metaTemp_",i, "_rr_0.txt", sep="")
            if (file.exists(fdir))
            {
                if (!file.info(fdir)$size == 0)
                {
                    Lt<-read.table(fdir,sep = "\t",header=F)
                    if (length(mL) == ncol(Lt))
                        names(Lt)<-mL 
                    if (nrow(mL) == ncol(Lt))
                        names(Lt)<-t(mL) 
                    L<-cbind(L, Lt)
                }
            }
            fdir <- paste(dirOP, "/beta_sam_",i, "_rr_0.txt", sep="")
            if (file.exists(fdir))
            {
                if (!file.info(fdir)$size == 0) 
                {
                    betasamInd<-read.table(fdir,sep = "\t", header=F)
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
                the<-cbind(the, read.table(fdir, sep = "\t", header=F))
            } 
            }
            fdir <- paste(dirOP,  "/metaFit_sam_",i, "_rr_0.txt", sep="")
            if (file.exists(fdir))
            {    
                if (!file.info(fdir)$size == 0) 
                {
                    sfitsam<-cbind(sfitsam, read.table(fdir, sep = "\t", header=F))
                }
            }
            fdir <- paste(dirOP, "/lambda_sam_",i,"_rr_0.txt", sep="")
            if (file.exists(fdir))
            {    
                if (!file.info(fdir)$size == 0) 
                {
                    lam<-cbind(lam,read.table(fdir,header=F))
                }  
            }
            fdir <- paste(dirOP, "/metaIndFit_sam_",i,"_rr_0.txt", sep="") 
            if (file.exists(fdir))
            {
                if (!file.info(fdir)$size == 0)
                {
                    metafitsamInd<-read.table(fdir,sep="\t",header=F) 
                    for (j in 1:bcol)
                        metaindfitsam<-cbind(metaindfitsam, t(t(metafitsamInd[,((j-1)*brow+1):(j*brow)])*betasamInd[,j]))                
                } 
            }
        }
    }
    ## column name for beta
    if (!is.null(bet))
    {
        if (ncol(bet) == ncol(specTitle))
            names(bet)<- t(specTitle[2,])
    }
    ## column name for delta
    if (!is.null(delmean))
    {
        if (ncol(delmean) == ncol(specTitle))
            names(delmean)<- t(specTitle[2,])
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
    
    fdir <- paste(dirOP, "/specFit_1_rr_1.txt", sep="")
	if (file.exists(fdir))
	{                
        if (!file.info(fdir)$size == 0) 
        {
            rrr<-read.table(fdir,sep = "\t",header=T)
        }
    }
    fdir <- paste(dirOP, "/NMRdata_mod_1.txt", sep="")
    if (file.exists(fdir))
    {                
        if (!file.info(fdir)$size == 0) 
        {
            Ndata<-read.table(fdir, header = FALSE, sep = "\t")
        }
    }
    fdir <- paste(dirOP, "/specFitHR_1_rr_",rr,".txt", sep="")
    if (file.exists(fdir))
    {                
        if (!file.info(fdir)$size == 0) 
        {
            rrrH<-read.table(fdir,sep = "\t",header=F)
            rrrH<-cbind(Ndata,rrrH)
        }
    }
    fdir <- paste(dirOP, "/metaTempHR_1_rr_",rr,".txt", sep="")
    if (file.exists(fdir))
    {                
        if (!file.info(fdir)$size == 0) 
        {
            LrrH<-read.table(fdir,sep = "\t",header=F)
            if (length(mL) == ncol(LrrH))
                names(LrrH)<-mL
            if ( nrow(mL) == ncol(LrrH))
                names(LrrH)<-t(mL)
        }
    }
    fdir <- paste(dirOP, "/metaTemp_1_rr_",rr,".txt", sep="")
    if (file.exists(fdir))
    {                
        if (!file.info(fdir)$size == 0) 
        {
            Lrr<-read.table(fdir,sep = "\t",header=F)
            if (length(mL) == ncol(Lrr))
               names(Lrr)<-mL
            if ( nrow(mL) == ncol(Lrr))
               names(Lrr)<-t(mL)
        }
    }
    fdir <- paste(dirOP, "/beta_1_rr_",rr,".txt", sep="")
    if (file.exists(fdir)) 
    {                
        if (!file.info(fdir)$size == 0) 
        {
            betrr<-read.table(fdir,header=F)
            if (length(mL) == nrow(betrr))
            row.names(betrr)<-mL
            if (nrow(mL) == ncol(betrr))
            row.names(betrr)<-t(mL)
        }
    }
    fdir <- paste(dirOP, "/beta_sam_1_rr_",rr,".txt", sep="")
    if (file.exists(fdir)) 
    {
        if (!file.info(fdir)$size == 0) 
        {
            betasamrr<-read.table(fdir,sep = "\t",header=F)
            if (length(mL) == nrow(betasamrr))
                row.names(betasamrr)<-mL
            if (nrow(mL) == nrow(betasamrr))
                row.names(betasamrr)<-t(mL)
        }
    }
    fdir <- paste(dirOP, "/theta_sam_1_rr_",rr,".txt", sep="")
    if (file.exists(fdir))
    {
        if (!file.info(fdir)$size == 0) 
        {
        	therr<-read.table(fdir,sep = "\t",header=F)
        }
    }
    fdir <- paste(dirOP, "/metaFit_sam_1_rr_",rr,".txt", sep="")
    if (file.exists(fdir))
    {if (!file.info(fdir)$size == 0) 
                {
      	sfitsamrr<-read.table(fdir,sep = "\t",header=F)
    }}
    ## individual metabolite fit posteriors
    brow<-nrow(betasamrr)
    bcol<-ncol(betasamrr)
    fdir <- paste(dirOP, "/metaIndFit_sam_1_rr_",rr,".txt", sep="")  
    if (file.exists(fdir))
    { 
        if (!file.info(fdir)$size == 0)
        {
             metafitsamrrInd<-read.table(fdir,sep="\t",header=F) 
             for (j in 1:bcol)
             {
                 if (j == 1)
                     metaindfitsamrr<-t(t(metafitsamrrInd[,((j-1)*brow+1):(j*brow)])*betasamrr[,j])
                 else
                     metaindfitsamrr<-cbind(metaindfitsamrr, t(t(metafitsamrrInd[,((j-1)*brow+1):(j*brow)])*betasamrr[,j]))
             }                       
        }
    }
    ## more spectra
    if (NoSpectra>1)
    {
        for (i in 2:NoSpectra)   
        { 
            fdir <- paste(dirOP, "/specFit_", i,"_rr_",rr,".txt", sep="")
            if (file.exists(fdir))
            {
                if (!file.info(fdir)$size == 0) 
                {
                rrr<-cbind(rrr, read.table(fdir,sep = "\t",header=T))
                }
            }           
            fdir <- paste(dirOP, "/NMRdata_mod_",i,".txt", sep="")
            if (file.exists(fdir))
            {
                if (!file.info(fdir)$size == 0) 
                {
                    Ndata<-read.table(fdir, header = FALSE, sep = "\t")
                }
            }
            fdir <- paste(dirOP, "/specFitHR_", i,"_rr_",rr,".txt", sep="")
            if (file.exists(fdir))
            {
                if (!file.info(fdir)$size == 0) 
                {
                    rrrH<-cbind(rrrH, Ndata, read.table(fdir,sep = "\t",header=F))
                }
            }
            fdir <- paste(dirOP, "/metaTempHR_",i, "_rr_",rr,".txt", sep="")
            if (file.exists(fdir))
            {
                if (!file.info(fdir)$size == 0) 
                {
                    LrrtH<-read.table(fdir,sep = "\t",header=F)
                    if (length(mL) == ncol(LrrtH))
                        names(LrrtH)<-mL
                    if ( nrow(mL) == ncol(LrrtH))
                        names(LrrtH)<-t(mL)
                    LrrH<-cbind(LrrH, LrrtH)
                    }
            }
            fdir <- paste(dirOP, "/metaTemp_",i, "_rr_",rr,".txt", sep="")
            if (file.exists(fdir))
            {
                if (!file.info(fdir)$size == 0) 
                {
                    Lrrt<-read.table(fdir,sep = "\t",header=F)
                    if (length(mL) == ncol(Lrrt))
                        names(Lrrt)<-mL
                    if ( nrow(mL) == ncol(Lrrt))
                        names(Lrrt)<-t(mL)
                    Lrr<-cbind(Lrr, Lrrt)
                }   
            }
            fdir <- paste(dirOP, "/beta_",i, "_rr_",rr,".txt", sep="")
            if (file.exists(fdir)) 
            {
                if (!file.info(fdir)$size == 0) 
                {
                    betrr<-cbind(betrr, read.table(fdir,header=F))
                }
            }            
            fdir <- paste(dirOP, "/beta_sam_",i, "_rr_",rr,".txt", sep="")
            if (file.exists(fdir)) 
            {
                if (!file.info(fdir)$size == 0) 
                {
                    betasamrrInd<-read.table(fdir,sep = "\t",header=F)
                    betasamrr<-cbind(betasamrr, betasamrrInd)  
                }
            }
            fdir <- paste(dirOP,  "/theta_sam_",i, "_rr_",rr,".txt", sep="")
            if (file.exists(fdir))
            {
                if (!file.info(fdir)$size == 0) 
                {
                    therr<-cbind(therr, read.table(fdir, sep = "\t", header=F))    
                }
            }
            fdir <- paste(dirOP,  "/metaFit_sam_",i, "_rr_",rr,".txt", sep="")
            if (file.exists(fdir))
            {
                if (!file.info(fdir)$size == 0) 
                {
                    sfitsamrr<-cbind(sfitsamrr, read.table(fdir, sep = "\t", header=F))    
                }   
            }       
            fdir <- paste(dirOP, "/metaIndFit_sam_",i, "_rr_",rr,".txt", sep="") 
            if (file.exists(fdir))
            {
                if (!file.info(fdir)$size == 0) 
                {
                    metafitsamrrInd<-read.table(fdir,sep="\t",header=F)   
                    for (j in 1:bcol)
                        metaindfitsamrr<-cbind(metaindfitsamrr, t(t(metafitsamrrInd[,((j-1)*brow+1):(j*brow)])*betasamrrInd[,j]))                    
                }
            }
        }    
    }
    ## column name for beta
    if (!is.null(betrr))
    {
        if (ncol(betrr) == ncol(specTitle))
            names(betrr)<- t(specTitle[2,])
		row.names(betrr)<-row.names(bet)

    }
    if (!is.null(bet))
    {
        bet2<-bet
        Metabolite<-row.names(bet2)
        row.names(bet2)<-NULL
        bet3<-cbind(Metabolite,bet2)
        write.table(bet3,file=paste(dirOP,"/RelCon.txt",sep=""),sep = "\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
    }
    if (!is.null(betrr))
    {
        betrr2<-betrr
        Metabolite<-row.names(betrr2)
        row.names(betrr2)<-NULL
        betrr3<-cbind(Metabolite,betrr2)
        write.table(betrr3,file=paste(dirOP,"/RelConRerun.txt",sep=""),sep = "\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
    }
    if (!is.null(delmean))
    {
        delmean2<-delmean
        Multiplet<-row.names(delmean2)
        row.names(delmean2)<-NULL
        delmean3<-cbind(Multiplet,delmean2)
        write.table(delmean3,file=paste(dirOP,"/MultipletsPpmShifts.txt",sep=""),sep = "\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
    }
    ns<-5
    betvA <-NULL  
    betvArr <-NULL 
    vA<-NULL
    vArr<-NULL
    if (!is.null(betasam))
    {
        f<-nrow(betasam)
        fc<-ncol(betasam)
        sno<-length(r)/ns
        ind<-fc/sno
        for (i in 1:sno) 
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
        write.table(vA2,file=paste(dirOP,"/RelConCreInt.txt",sep=""),sep = "\t",row.names = TRUE,col.names = FALSE,quote=FALSE)
    }
    if (!is.null(betasamrr))
    { 
        f<-nrow(betasamrr)
        fc<-ncol(betasamrr)
        ind<-fc/sno
        for (i in 1:sno) 
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
        write.table(vArr2,file=paste(dirOP,"/RelConCreIntRerun.txt",sep=""),sep = "\t",row.names = TRUE,col.names = FALSE,quote=FALSE)
    }   
	## batman fitting results
    BM<-list(specTitle = specTitle, sFit=r,sFitHR=rH, beta=bet, betaSam=betasam, betaCI=vA, metaTemp=L,
    metaTempHR=LH, metaFitSam = sfitsam, metaIndFitSam = metaindfitsam, thetaSam=the, delta = delmean, deltaSam = del,
    sFitRerun=rrr, sFitRerunHR=rrrH, betaRerun = betrr, betaSamRerun = betasamrr, betaCIRerun=vArr,
    metaTempRerun = Lrr, metaTempRerunHR = LrrH, metaFitSamRerun=sfitsamrr, metaIndFitSamRerun=metaindfitsamrr,
    thetaSamRerun=therr,outputDir = dirOP) 
    warnRead<-options(warn = warnDef)
    return (BM)
}    
                   


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
	    dirA<-c(extdir,extdir,extdir)
	}
	ctime <- format(Sys.time(), "%d_%b_%H_%M_%S")

	dir1<-paste(dirA[2],"/batmanOptions.txt",sep="")
	dir2<-paste(dirA[2],"/NMRdata.txt",sep="")
	dir3<-paste(dirA[2],"/metabolitesList.txt",sep="")
	dir4<-paste(dirA[2],"/multi_data.dat",sep="")
	
    dirctime<-paste(dirA[3],"/",ctime,sep="")
	if(!file.exists(dirctime)) {
    	dir.create(dirctime)
	}

	dir5<-paste(dirctime,"/",sep="")
	
	cp<-file.copy(dir1,dir5)
	dirList<-paste(dirA[2],"/metabolitesList.csv",sep="")
	cp<-file.copy(dirList,dir5)
	
	filedir<-c(dir1,dir2,dir3,dir4,dir5)
	
	dirL<-paste(dirA[2],"/metabolitesList.csv",sep="")
	dirR<-paste(dirA[2],"/multi_data.csv",sep="")
	dirRU<-paste(dirA[2],"/multi_data_user.csv",sep="")
    
    ## read in batman optitons
    con  <- file(dir1, open = "r")
    oneLine <- readLines(con, n = 30, warn = FALSE)
    fL<-substr(oneLine,1,1)
    nL<-which(is.na(match(fL,"%")))
    myVector <- strsplit(oneLine[nL[2]], ":")
    sno <- as.numeric(myVector[[1]][2])
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
   
    if ((ncol(sa)-1)<sno)
        return(cat("No. of spectra included smaller than input spectra.\n"))
    if (!is.null(colnames(sa))) {
        saname<-colnames(sa)
        saname<-saname[2:length(saname)]
        stit<-rbind(1:length(saname),saname)
    } else {
        stit<-rbind(1:(ncol(sa)-1),1:(ncol(sa)-1))
    }
    ## spectra number
    write.table(stit[,1:sno],file=paste(dirctime, "/spectraTitle.txt", sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep = "\t")
    ## choose whether to parallelize spectra
    if (sno>1 && fixeff == 0)    
    {
		cat("\nHow many parallel processes (multicores) do you want to run the multi-spectra analysis?")
		cat("\n(Enter 1 for running them sequentially.)\n")
		cat("\n Parallel processing of multi spectra currently cannot display the progress\n")
        cat(" bar (or any words), if you input is > 1, please be patient for the results :)\n\n")
		wr<- getinput(lowlim=1,highlim=20)  
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
        	multispec<- foreach (nospec = 1:sno,.packages="batman") %dopar% {
                pBar <- txtProgressBar(min =0, max = ito, style = 3)
                out<-.Call("batman",filedir,as.integer(bn),as.integer(ito),as.integer(rr),as.integer(nospec-1), pBar, PACKAGE = "batman")
                close( pBar )
	}})[3]} else if (wr == 1 && fixeff == 1) { 
        stime <- system.time ({
        	multispec<- for (nospec in 1:1) {
                pBar <- txtProgressBar(min =0, max = ito, style = 3)
            	out<-.Call("batman",filedir,as.integer(bn),as.integer(ito),as.integer(rr),as.integer(sno-1),pBar,PACKAGE = "batman")
            	close( pBar )
	}})[3]} else { 
        stime <- system.time ({
        	multispec<- for (nospec in 1:sno) {
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
	BM<-readBatmanOutput(dirctime) 
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
	cat("\nCompleted.\n")
	warnRead<-options(warn = warnDef)
	return(BM)
}
