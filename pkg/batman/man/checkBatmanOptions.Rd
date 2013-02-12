\name{checkBatmanOptions}
\alias{checkBatmanOptions}
\title{Check batmanOptions.txt file for version later than 1.0.4}
\description{Check batmanOptions.txt file and add a new line input at the end of the file.}
\usage{checkBatmanOptions(dir)}
\arguments{
\item{dir}{The directory of batmanOptions.txt file.}
}
\examples{
library(batman)
## createfolder "runBATMAN" in current working directory
batmanDir = newDir(runBATMANDir = getwd(), overwriteFile = TRUE)
checkBatmanOptions(dir = paste(batmanDir[2], "/batmanOptions.txt", sep = ""))
}
\keyword{datasets}
