\name{createChemShiftPerSpec}
\alias{createChemShiftPerSpec}
\title{Creating the file chemShiftPerSpec.csv which contains chemical shift parameters for all multiplets and spectra.}
\description{This function creates a file called chemShiftPerSpec.csv, so user can specify chemical shift parameter for each spectrum and multiplet.  
 The first column is multiplet names in the same order as the template inputs in multi_data.csv and/or multi_data_user.csv 
 (depending on user choice of using one or both of them) file(s). The second column is the default chemical shift value (pos_in_ppm) for the corresponding multiplet.
 From the third column forward is the chemical shift value for each spectrum in the same order as they read in by BATMAN, if 'n' is present in the field,
 the default chemical shift value (or overwrite_pos value if given) will be used. }
\usage{
 createChemShiftPerSpec(tempOption, dirIP)}
\arguments{
\item{tempOption}{Choose template file(s). \code{tempOption = 1} for multi_data.csv, \code{tempOption = 2} for multi_data_user.csv, and
 \code{tempOption = 3} for both files.}
\item{dirIP}{The input directory of BATMAN. This is the path ending with '/BatmanInput' if runBATMAN directory is created.}
}
\seealso{
\code{\link{batman}}
}
\examples{
library(batman)
## createfolder "runBATMAN" in current working directory
batmanDir = newDir(runBATMANDir = getwd(), overwriteFile = TRUE) 
## create chemShiftPerSpec.csv
createChemShiftPerSpec(tempOption = 1, dirIP = batmanDir[2]) 
}
\keyword{aplot}
