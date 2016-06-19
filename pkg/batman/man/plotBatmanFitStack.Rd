\name{plotBatmanFitStack}
\alias{plotBatmanFitStack}
\title{Stack plot Batman Metabolite Fit of NMR Spectra (With Down Sampling)}
\description{This function plots the BATMAN fit results in stack, and saves the 
 figure to pdf file in specified directory. User can choose to plot all or some 
 of the spectra analyzed, show the metabolite fit and give a range for x and y limits.
 The figure file will not be overwritten if it already exists by default. A prefix can be added to 
 the file name for new saves.}
\usage{
plotBatmanFitStack(BM, offset = 1, mirroredWav = TRUE, specNo,  xfrom,  
                   xto, yfrom, yto, listMeta = FALSE, metaName,  
                   saveFig = TRUE, saveFigDir = BM$outputDir,  
                   prefixFig, rerun = FALSE, placeLegend = "topright", 
                   plotColour, overwriteFig = FALSE, 
                   metaLwd = 2, metaLty = 5, orientation = "L",
                   showPlot)
}
\arguments{
\item{BM}{\code{batman} output data frame.}
\item{offset}{Offset value for the stack plot.}
\item{mirroredWav}{Plot mirrored wavelet fit if set TRUE. The default is TRUE}
\item{specNo}{Vector of spectra ID in the input order to specify which spectra to be plotted in stack plot.}
\item{xfrom}{The start ppm value to plot. Default is set to the start ppm value of the whole processed range.}
\item{xto}{The end ppm value to plot. Default is set to the end ppm value of the whole processed range.}
\item{yfrom}{The start value of vertical axis to plot. Default is set to 0.}
\item{yto}{The end value of vertical axis to plot. Default is set to the maximum value of the spectrum point in display.}
\item{listMeta}{Individual metabolite fit will also be shown in the plot if set TRUE. The default is FALSE.}
\item{metaName}{One or more specified metabolite fits will be shown in the plot.
 If no name was given and \code{listMeta = TRUE}, all the individual metabolite fit will be shown.}
\item{saveFig}{Save figure to pdf file if set TRUE. The default is TRUE.}
\item{saveFigDir}{Save figure in this directory. The default is output directory of BM.}
\item{prefixFig}{Add prefix to each saved figure name. The default is no prefix.}
\item{rerun}{Set to FALSE to plot \code{batman} result, and TRUE to plot \code{batmanrerun} result.}
\item{placeLegend}{Where to place the legend in figure. The default is "topright".}
\item{plotColour}{User can specify colours for each metabolite if \code{listMeta = TRUE}. If not, a set of randomly generated colours will be used.}
\item{overwriteFig}{Overwrite saved figure file in pdf format if \code{overwriteFig = TRUE}. If set to \code{FALSE}, a new figure file with system time as postfix will be created. The default is \code{FALSE}.}
\item{metaLwd}{The line widths for metabolite fit.}
\item{metaLty}{The line types for metabolite fit.}
\item{orientation}{The orientation of plot, either portrait or landscape, the defult is "L".}
\item{showPlot}{If set \code{FALSE}, no plot will be shown on display, the pdf file(s) for 
 the figure plot(s) will be created in output folder. If missing from input, for windows and osx operating systems, it will be set to \code{TRUE}, for the rest operating systems, it will be set to \code{FALSE} automatically.}
}
\examples{
library(batman)
## Run BATMAN 
if(interactive()) 
{
bm<-batman()
## then plot results
plotBatmanFitStack(bm)
}
}
\keyword{aplot}
