\name{readBruker}
\alias{readBruker}
\title{Read Binary Bruker NMR Spectra}
\description{
 Read in multiple binary Bruker NMR spectra (1D) of the same ppm scale from a specified 
 folder, and return a matrix with columns: \deqn{[ppm, spectrum 1, spectrum 2, ...].}}
\usage{readBruker(BrukerDataDir)}
\arguments{
\item{BrukerDataDir}{The directory of the folder containing 1D Bruker spectral data files. 
 Recursively finds all the "1r" files in \code{datapath} and read in.}
}
\value{It returns a matrix with columns: \deqn{[ppm, spectrum 1, spectrum 2, ...].}}
\examples{
library(batman)
## Read in all burker NMR spectra files, replace "your\data\path\here" with the
## directory of the data files you want to read.
## brukerdata<-readBruker("your\data\path\here")
}
\keyword{datasets}

 
