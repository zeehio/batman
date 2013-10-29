\name{plotShift}
\alias{plotShift}
\title{Boxplot or Histogram of ppm Shift Posterior distributions for Multiplets of Named Metabolite}
\description{This function provides boxplots or histograms of the ppm shift posterior distributions of multiplets, 
 and saves the figure to pdf file in specified directory. The file name is in the format
 of "spec_\eqn{i}_\eqn{metaName}_ppmShift.pdf", where \eqn{i} is the spectrum number and "metaName" is 
 the input metabolite name if given. The figure file will not be overwritten if 
 it already exists. A prefix can be given to the file name for new saves.}
\usage{
 plotShift(BM, metaName, plotHist = FALSE, breaks, perMult = FALSE,
           saveFig = TRUE, saveFigDir = BM$outputDir, prefixFig, 
           overwriteFig = FALSE)}
\arguments{
\item{BM}{\code{batman} output data frame.}
\item{metaName}{Only multiplets belonging to the named Metabolite will be shown. 
 Only one metabolite name can be given. If missing, all metabolites will be plotted.}
\item{plotHist}{If \code{plotHist = TRUE}, the ppm shift posteriors will be displayed as histogram. 
 The default is FALSE.}
\item{breaks}{A single number to set the number of bins for the histogram. If missing from the input, 
 it is set to the data length divided by 3.} 
\item{perMult}{If set TRUE plot the shifts per multiplet, otherwise, plot the shifts per spectrum.}
\item{saveFig}{Save figure to pdf file if set TRUE. The default is TRUE.}
\item{saveFigDir}{Save pdf file in this directory. The default is the output directory of BM.}
\item{prefixFig}{Add prefix to each saved figure name. The default is no prefix.}
\item{overwriteFig}{Overwrite saved figure file in pdf format if \code{overwriteFig = TRUE}. The default is \code{FALSE}}.
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
## Plot ppm shift for each multiplet.
plotShift(bm)
}
}
\keyword{aplot}
