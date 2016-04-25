\name{plotBatmanFitHR}
\alias{plotBatmanFitHR}
\title{Plot BATMAN Metabolite Fit of NMR Spectra in Original Resolution (Without Down Sampling)}
\description{This function plots a high resolution BATMAN fit results 
 (without down sampling), and save figure 
 to pdf file in user specified directory. For multiple spectra analysis, 
 the file name is in the format 
 of "specFitHR_\eqn{i}_\eqn{metaName}.pdf", where \eqn{i} is the spectrum number 
 in the figure and the metabolite name 
 will be shown in place of \eqn{metaName} if supplied. 
 The figure file will not be overwritten if it already exists. 
 A prefix can be given to the file name for new saves.}
\usage{
 plotBatmanFitHR(BM, xfrom, xto, yfrom, yto, metaName, saveFig = TRUE, 
                 saveFigDir = BM$outputDir, prefixFig, rerun = FALSE, 
                 overwriteFig = FALSE, showPlot = TRUE)}
\arguments{
\item{BM}{\code{batman} output data frame.}
\item{xfrom}{The start ppm value to plot. Default is set to the start ppm value of the whole processed range.}
\item{xto}{The end ppm value to plot. Default is set to the end ppm value of the whole processed range.}
\item{yfrom}{The start value of vertical axis to plot. Default is set to 0.}
\item{yto}{The end value of vertical axis to plot. Default is set to the maximum value of the spectrum point in display.}
\item{metaName}{Individual metabolite fit will also be shown in the plot if a metabolite name is given.
 Only one metabolite name can be given, if missing from input all metabolites will be plotted.}
\item{saveFig}{Save figure(s) to pdf file(s) if set TRUE. The default is TRUE.}
\item{saveFigDir}{Save figure(s) in this directory. The default is the output directory of BM.}
\item{prefixFig}{Add prefix to each saved figure name. The default is no prefix.}
\item{rerun}{Set to FALSE to plot \code{batman} result, and TRUE to plot \code{batmanrerun} result.}
\item{overwriteFig}{Overwrite saved figure file in pdf format if \code{overwriteFig = TRUE}. The default is \code{FALSE}.}
\item{showPlot}{If set \code{FALSE}, no plot will be shown on display, the pdf file(s) for 
 the plot(s) will be created in output folder. The default is \code{TRUE}.}
}
\seealso{ 
\code{\link{batman}, \link{batmanrerun}}
}
\examples{
library(batman)
## Run BATMAN fit
if(interactive()) 
{
bm<-batman()
## Plot batman Fit in its original resolution if the option parameter  
## is set to 1 for "Save metabolites fit same as the original spectrum  
## resolution (1/0)" in "batmanOptions.txt", .
plotBatmanFitHR(bm)
}
}
\keyword{aplot}
