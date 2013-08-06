\name{plotDiagnosticScatter}
\alias{plotDiagnosticScatter}
\title{Diagnostic scatter plot of batman metabolites fit vs NMR spectra bins or minimum wavelet fit.}
\description{When fitting a large number of spectra, this plot facilitates discovery of spectra or 
 metabolites which are poorly fit.}
\usage{
plotDiagnosticScatter(BM, binWidth = 0.018, cexID = 0.5, saveFig = TRUE, 
                      saveFigDir = BM$outputDir, prefixFig, 
                      rerun = FALSE, placeLegend = "topright", 
                      overwriteFig = FALSE)
}
\arguments{
\item{BM}{\code{batman} output data frame.}
\item{binWidth}{The full width of the bins to integrate. The centre of a bin is the estimated mean posterior 
 chemical shift for each multiplet in each spectrum.}
\item{cexID}{Character size for the spectra ID number.}
\item{saveFig}{Save figure to pdf file if set TRUE. The default is TRUE.}
\item{saveFigDir}{Save figure in this directory. The default is output directory of BM.}
\item{prefixFig}{Add prefix to each saved figure name. The default is no prefix.}
\item{rerun}{Set to FALSE to plot \code{batman} result, and TRUE to plot \code{batmanrerun} result.}
\item{placeLegend}{Where to place the legend in figure. The default is "topright".}
\item{overwriteFig}{Overwrite saved figure file in pdf format if \code{overwriteFig = TRUE}. The default is \code{FALSE}}
}
\examples{
library(batman)
## Run BATMAN 
if(interactive()) 
{
bm<-batman()
## then plot results
plotDiagnosticScatter(bm)
}
}
\keyword{aplot}
