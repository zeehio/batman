\name{checkBatmanOptions}
\alias{checkBatmanOptions}
\title{Check previouse versions of batmanOptions.txt file and unify the parameter
names to the current one.}
\description{Check batmanOptions.txt file and may add a new input line at the 
end of the file for old versions.}
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
