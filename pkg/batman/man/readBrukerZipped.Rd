\name{readBrukerZipped}
\alias{readBrukerZipped}
\title{Read Raw Binary Bruker NMR Spectra in Zipped format}
\description{
 Read in multiple raw binary Bruker NMR spectra (1D), with spectrum data in a zipped format, 
 from a specified folder, and return a matrix with columns: \deqn{[ppm, spectrum 1, spectrum 2, ...].}
 Interpolation may be performed if spectra have different ppm scales.}
\usage{readBrukerZipped(BrukerDataZipDir)}
\arguments{
\item{BrukerDataZipDir}{The directory of the folder containing zipped 1D Bruker spectral data files. 
 Recursively finds all the "*.zip" files in \code{datapath}, unzipped them in the same folder, call
 "*.zip" files to read in spectra, and delete the unipped folders. If no "*.zip" file was found,
 it works the same as "*.zip" files.}
}
\value{It returns a matrix same as \code{readBruker}, 
 with columns: \deqn{[ppm, spectrum 1, spectrum 2, ...].}}
\examples{
library(batman)
## Read in all Burker NMR spectra files, replace "/your/data/path/here" with the
## directory of the data files you want to read.
## brukerdata<-readBrukerZipped("/your/data/path/here")
}
\keyword{datasets}

 
