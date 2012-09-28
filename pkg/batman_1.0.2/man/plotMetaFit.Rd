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
             saveFigDir = BM$outputDir, prefixFig, rerun = FALSE)	
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
}
\seealso{
 \code{\link{batman}, \link{batmanrerun}}
}
\examples{
library(batman)
## Run BATMAN fit, then plot metabolite fit
bm<-batman()
## Plot metabolites Fit.
plotMetaFit(bm)
}
\keyword{aplot}
