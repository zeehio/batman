\name{plotRelCon}
\alias{plotRelCon}
\title{Boxplot or Histogram of Posterior distributions of Relative Concentrations 
 for Listed Metabolites with 95\% Credible Interval}
\description{This function plots the posterior distributions of relative concentrations, and saves the figure 
 to pdf file. The file name is in the format of "spec_\eqn{i}_RelCon_\eqn{j1}to\eqn{j2}.pdf", 
 where \eqn{i} are the spectrum numbers and \eqn{j1} and \eqn{j2} are the order 
 numbers of fitted metabolites in the order of their input in file metaboliteList.csv. 
 The figure file will not be overwritten if it already exists. 
 A prefix can be added to file name for new saves.}
\usage{
 plotRelCon(BM, metaName, plotHist = FALSE, breaks, 
            saveFig = TRUE, saveFigDir = BM$outputDir, 
            prefixFig, rerun = FALSE, overwriteFig = FALSE)}
\arguments{
\item{BM}{\code{batman} output data frame.}
\item{metaName}{Only multiplets belonging to the named Metabolite will be shown. 
 Only one metabolite name can be given. If missing, all metabolites will be plotted.}
\item{plotHist}{If \code{plotHist = TRUE}, the ppm shift posteriors will be displayed as histogram. 
 The default is FALSE.}
\item{breaks}{A single number to set the number of bins for the histogram. If missing from the input, 
 it is set to the data length divided by 3.} 
\item{saveFig}{Save figure(s) to pdf file(s) if set TRUE. The default is TRUE.}
\item{saveFigDir}{Save figure(s) in this directory. The default is output directory of BM.}
\item{prefixFig}{Add prefix to each saved figure name. The default is no prefix.}
\item{rerun}{Set to FALSE to plot \code{batman} result, and TRUE to plot \code{batmanrerun} result.}
\item{overwriteFig}{Overwrithe saved figure file in pdf format if \code{overwriteFig = TRUE}. The default is \code{FALSE}}.
}
\seealso{
\code{\link{batman}, \link{batmanrerun}}
}
\examples{
library(batman)
## Run BATMAN and then plot relative concentration
bm<-batman()
100
1
## Plot relative concentrations
plotRelCon(bm)
}
\keyword{aplot}
