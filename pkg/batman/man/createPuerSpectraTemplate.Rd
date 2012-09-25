\name{createPureSpectraTemplate}
\alias{createPureSpectraTemplate}
\title{Creating a folder called 'PureSpectraTemplate' in the specified input directory. The folder contains pure metabolite spectrum template in .txt file with 
metabolite name as the file name.}
\description{This function will read in pure metabolites spectra in Bruker format and save them in .txt format in folder PureSpectraTemplate. This folder will be used if the 'couple_code' value in multi_data.csv and/or multi_data_user.csv is set to '-2'.}
\usage{
 createPureSpectraTemplate(dirPureSpec, metaNames, dirIP)}
\arguments{
\item{dirPureSpec}{A vector containing the directories of Bruker pure metabolite spectra files.}
\item{metaNames}{The vector of metabolites names in the same order as the spectra directories in \code{dirPureSpec}.}
\item{dirIP}{The input directory of BATMAN. This is the path ending with '/runBATMAN/BatmanInput' if runBATMAN directory is created.}
}
\examples{
library(batman)
## createfolder "runBATMAN" in current working directory
batmanDir = newDir(runBATMANDir = getwd(), overwriteFile = TRUE) 
## create chemShiftPerSpec.csv
createPureSpectraTemplate(dirPureSpec = c(paste(batmanDir[2],"/testBrukerPureSpec/testPure1",sep = ""),
paste(batmanDir[2],"/testBrukerPureSpec/testPure2",sep = "")), metaNames = c("testPure1", "testPure2"), 
dirIP = batmanDir[2])
}
\keyword{aplot}
