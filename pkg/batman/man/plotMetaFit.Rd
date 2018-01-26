\name{plotMetaFit}
\alias{plotMetaFit}
\title{Plot Posterior Means of Metabolites Fit with 95\% Credible Interval}
\description{This function plots posterior means of the metabolite fit with 95\% 
 credible interval , and saves the figure 
 to pdf file in specified directory. For multiple metabolites, the file name is in 
 the format of "spec_\eqn{i}to\eqn{j}_mFitSam.pdf", where \eqn{i} and \eqn{j} 
 are range numbers of spectra in the figure. 
 A maximum of 2 spectra will be shown in each figure. Figure file 
 will not be overwritten if it already exists.  Prefix can be added to the file name for new saves.}
\usage{
 plotMetaFit(BM, from, to, metaName, saveFig = TRUE,
             saveFigDir = BM$outputDir, prefixFig, 
             rerun = FALSE, overwriteFig = FALSE,
             showPlot)	
}
\arguments{
\item{BM}{\code{batman} output data frame.}
\item{from}{The start ppm value to plot. Default is set to the start ppm value of the whole processed range.}
\item{to}{The end ppm value to plot. Default is set to the end ppm value of the whole processed range.}
\item{metaName}{Only multiplets belonging to the named Metabolite will be shown. 
 Only one metabolite name can be given. If missing, all metabolites will be plotted.}
\item{saveFig}{Save figure to pdf file if set TRUE. The default is TRUE.}
\item{saveFigDir}{Save figure in this directory. The default is current working directory.}
\item{prefixFig}{Add prefix to each saved figure name. The default is no prefix.}
\item{rerun}{Set to FALSE to plot \code{batman} result, and TRUE to plot \code{batmanrerun} result.}
\item{overwriteFig}{Overwrite saved figure file in pdf format if \code{overwriteFig = TRUE}. If set to \code{FALSE}, a new figure file with system time as postfix will be created. The default is \code{FALSE}.}
\item{showPlot}{If set \code{FALSE}, no plot will be shown on display, the pdf file(s) for 
 the figure plot(s) will be created in output folder. If missing from input, for windows and osx operating systems, it will be set to \code{TRUE}, for the rest operating systems, it will be set to \code{FALSE} automatically.}
}
\seealso{
 \code{\link{batman}, \link{batmanrerun}}
}
\examples{
library(batman)
## Run BATMAN fit, then plot metabolite fit
if(interactive()) 
{
bm<-batman()
## Plot metabolites Fit.
plotMetaFit(bm)
}
}
\keyword{aplot}
