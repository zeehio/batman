\name{readBatmanOutput}
\alias{readBatmanOutput}
\title{Reads in BATMAN Output Data Files}
\description{Reads in output data files from \code{batman} in specified folder.}
\usage{readBatmanOutput(dirOP,dirIP)}
\arguments{
\item{dirOP}{The folder with \code{batman} output files.}
\item{dirIP}{The folder with \code{batman} input files.}
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
bm<-batman()
100
1
## Read in output files in saved directory.
bmread<-readBatmanOutput(bm$outputDir,bm$inputDir)
}
\keyword{datasets}
