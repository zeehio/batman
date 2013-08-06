\name{plotChemShiftDist}
\alias{plotChemShiftDist}
\title{Plot histogram of chemical shifts for the multiplets across a series of spectra.}
\description{This function plots the histogram of the mean posterior estimated chemical 
 shifts for the multiplets of certain or all metabolites across a series of spectra. User can choose 
 to plot all or some of the metabolite. The figure file will not be overwritten if it 
 already exists by default. A prefix can be added to the file name for new saves.}
\usage{
plotChemShiftDist(BM, metaName, breaks = 20, xlim, 
                  saveFig = TRUE, saveFigDir = BM$outputDir,
                  prefixFig, overwriteFig = FALSE)
}
\arguments{
\item{BM}{\code{batman} output data frame.}
\item{metaName}{One or more specified metabolites will be shown.
 If no name was given, all the individual metabolites will be shown.}
\item{breaks}{A single number to set the number of bins for the histogram.} 
\item{xlim}{The range of x values.}
\item{saveFig}{Save figure to pdf file if set TRUE. The default is TRUE.}
\item{saveFigDir}{Save figure in this directory. The default is output directory of BM.}
\item{prefixFig}{Add prefix to each saved figure name. The default is no prefix.}
\item{overwriteFig}{Overwrite saved figure file in pdf format if \code{overwriteFig = TRUE}. The default is \code{FALSE}}
}
\examples{
library(batman)
## Run BATMAN 
if(interactive()) 
{
bm<-batman()
## then plot results
plotChemShiftDist(bm)
}
}
\keyword{aplot}
