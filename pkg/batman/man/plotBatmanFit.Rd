\name{plotBatmanFit}
\alias{plotBatmanFit}
\title{Plot Batman Metabolite Fit of NMR Spectra (With Down Sampling)}
\description{This function plots the BATMAN fit results, and saves the 
 figure to pdf file in specified directory. For multiple spectra analysis, 
 the file name is in the format of "specFit_\eqn{i}to\eqn{j}_\eqn{metaName}.pdf", where \eqn{i} and \eqn{j} are 
 the range numbers of spectra in the figure and the metabolite name will be shown 
 in place of \eqn{metaName} if supplied. Maximum of 2 spectra will be shown in each figure. 
 The figure file will not be overwritten if it already exists by default. A prefix can be added to 
 the file name for new saves.}
\usage{
 plotBatmanFit(BM, xfrom, xto, yfrom, yto, listMeta = FALSE,
               metaName, saveFig = TRUE, saveFigDir = BM$outputDir, 
               prefixFig, rerun = FALSE, placeLegend,
               plotColour, overwriteFig = FALSE,
               showPlot = TRUE)}
\arguments{
\item{BM}{\code{batman} output data frame.}
\item{xfrom}{The start ppm value to plot. Default is set to the start ppm value of the whole processed range.}
\item{xto}{The end ppm value to plot. Default is set to the end ppm value of the whole processed range.}
\item{yfrom}{The start value of vertical axis to plot. Default is set to 0.}
\item{yto}{The end value of vertical axis to plot. Default is set to the maximum value of the spectrum point in display.}
\item{listMeta}{Individual metabolite fit will also be shown in the plot if set TRUE. The default is FALSE.}
\item{metaName}{One or more specified metabolite fits will be shown in the plot.
 If no name was given and \code{listMeta = TRUE}, all the individual metabolite fit will be shown.}
\item{saveFig}{Save figure(s) to pdf file(s) if set TRUE. The default is TRUE.}
\item{saveFigDir}{Save figure(s) in this directory. The default is output directory of BM.}
\item{prefixFig}{Add prefix to each saved figure name. The default is no prefix.}
\item{rerun}{Set to \code{FALSE} to plot \code{batman} result, and TRUE to plot \code{batmanrerun} result.}
\item{placeLegend}{Where to place the legend in figure. The default is "topright".}
\item{plotColour}{User can specify colours for each metabolite if \code{listMeta = TRUE}. If not, a set of randomly generated colours will be used.}
\item{overwriteFig}{Overwrite saved figure file in pdf format if \code{overwriteFig = TRUE}. The default is \code{FALSE}.}
\item{showPlot}{If set \code{FALSE}, no plot will be shown on display, the pdf file(s) for 
 the plot(s) will be created in output folder. The default is \code{TRUE}.}
}
\seealso{
\code{\link{batman}, \link{batmanrerun}}
}
\examples{
library(batman)
## Run BATMAN 
if(interactive()) 
{
bm<-batman()
## then plot results
plotBatmanFit(bm)
}
}
\keyword{aplot}
