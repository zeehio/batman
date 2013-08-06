\name{readBatmanOutput}
\alias{readBatmanOutput}
\title{Reads in BATMAN Output Data Files}
\description{Reads in output data files from \code{batman} in specified folder.}
\usage{readBatmanOutput(dirOP, dirIP,readMetaIndFitSam = TRUE,
                           readMetaTempHR = TRUE, readMetaTemp = TRUE)}
\arguments{
\item{dirOP}{The folder with \code{batman} output files.}
\item{dirIP}{The folder with \code{batman} input files.}
\item{readMetaIndFitSam}{If set \code{TRUE}, read in the posterior samples of individual metabolites fit.}
\item{readMetaTempHR}{If set \code{TRUE}, read in the posterior means of fitted metabolite templates 
 (without down sample).}
\item{readMetaTemp}{If set \code{TRUE}, read in the posterior means of fitted metabolite templates (down sampled).}
}
\value{	
 It returns a data list with the objects described in \code{\link{batman}}.
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
## Read in output files in saved directory.
bmread<-readBatmanOutput(bm$outputDir,bm$inputDir)
}
}
\keyword{datasets}
