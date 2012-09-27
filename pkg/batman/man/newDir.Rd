\name{newDir}
\alias{newDir}
\title{Create a New Directory with Installed Data Files}
\description{A folder name "runBATMAN" will be created in specified directory,
 within which, two folders "BatmanInput" and "BatmanOutput" are created.
 "BatmanInput" contains the input data files from installed package folder
 "extdata" and user inputted data files. The user only needs to modify files in this
 folder to change the settings for running \code{batman}. The \code{batman} 
 output files will be saved in "BatmanOutput" subfolders.}
\usage{newDir(runBATMANDir = getwd(), overwriteFile = FALSE)}
\arguments{
\item{runBATMANDir}{User specified directory for working with BATMAN, the default is current work directory. }
\item{overwriteFile}{If folder "runBATMAN" exist, set TRUE will overwrite folder. 
 The default is FALSE.}
}
\seealso{	
\code{\link{batman}}
}
\examples{
library(batman)
## This will create the folder "runBATMAN" in current working directory, 
## within the folder "runBATMAN", a subfolder "BatmanInput" contains all the
## input files batman uses. Users can modify the input files to change 
## the settings for running batman.
newDir(runBATMANDir = getwd(), overwriteFile = TRUE) 
}
\keyword{internal}
