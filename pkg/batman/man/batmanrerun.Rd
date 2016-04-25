\name{batmanrerun}
\alias{batmanrerun}
\title{Perform BATMAN with Fixed (Previously Estimated) Multiplet Positions}
\description{This performs metabolite and wavelet fitting to input NMR spectra with
 fixed multiplet position obtained from running \code{batman}, and also 
 plots fitting results. The user should modify parameters in the copy file "batmanOptions.txt" 
 in \code{batman} output folder to change the rerun settings. }
\usage{
 batmanrerun(BM, figBatmanFit = TRUE, listMeta = FALSE, 
             figRelCon = FALSE, figMetaFit = FALSE, showPlot = TRUE)
 }
\arguments{
\item{BM}{\code{batman} output data frame.}
\item{figBatmanFit}{Plot metabolites and wavelets fit if set \code{TRUE}. The default is \code{TRUE}.}
\item{listMeta}{Individual metabolite fit will also be shown in the plot if set \code{TRUE}. The default is \code{FALSE}.}
\item{figRelCon}{Plot posterior samples of the relative concentration for listed 
 metabolites with 95\% credible interval if set \code{TRUE}. The default is \code{FALSE}.}
\item{figMetaFit}{If set \code{TRUE}, plot the posterior mean of the metabolites fit with 95\% 
 credible interval. The default is \code{FALSE}.}
\item{showPlot}{If set \code{FALSE}, no plot will be shown on display, the pdf file(s) for 
 the plot(s) will be created in output folder. The default is \code{TRUE}.}
}
\value{
 When \code{batmanrerun} is called with multiplet ppm shifts fixed from the 
 \code{batman} results, the following objects are added to the \code{batman} result:
 
 \item{sFitRerun}{A matrix \eqn{t \times 5n} of BATMAN rerun fit results (down sampled). 
 For 1 spectrum, it is a matrix with 5 columns: \deqn{[ppm, original spectrum, 
 metabolites fit, wavelets fit, overall fit].} \eqn{n} is the number of spectra,
 and \eqn{t} is the number of data points in each spectrum.} 
 
 \item{sFitRerunHR}{A matrix \eqn{t \times 3n} of BATMAN rerun fit results in the 
 original resolution (without down sample). For 1 spectrum, it is a matrix with 3 columns:  
 \deqn{[ppm, original spectrum, metabolites fit].} \eqn{n} is the number of spectra, and
 \eqn{t} is the number of data points (without down sample) in each spectrum.} 
 
 \item{betaRerun}{For batman rerun, a matrix (\eqn{m \times n}) containing the 
 posterior means of relative concentrations for \eqn{m} fitted metabolites and 
 \eqn{n} spectra.}
 
 \item{betaSamRerun}{For batman rerun, a matrix (\eqn{m \times (s*n)}) containing 
 (for the first spectrum) \eqn{s} posterior samples of the relative concentrations 
 in its rows. \eqn{m} is the number of fitted metabolites. \eqn{n} is the number of spectra analyzed. 
 The subsequent columns contain the same data format for the rest \eqn{n-1} spectra.}
  
 \item{betaCIRerun}{For batman rerun, a matrix (\eqn{m \times 2n}) containing the
 95\% credible interval of the relative concentrations for \eqn{m} fitted metabolites. 
 Every pair of columns is for one spectrum.}
 
 \item{metaTempRerun}{For batman rerun, a matrix (\eqn{t \times (m*n)}) containing 
 the posterior means of \eqn{m} fitted metabolite templates in its columns (down sampled).
 \eqn{n} is the number of spectra analyzed and \eqn{t} is the number of data 
 points in each spectrum.}
 
 \item{metaTempRerunHR}{For batman rerun, a matrix (\eqn{t \times (m*n)}) 
 containing the posterior means of \eqn{m} fitted metabolite templates in its 
 columns (without down sample). \eqn{n} is the number of spectra analyzed and 
 \eqn{t} is the number of data points (without down sample) in each spectrum.}
     
 \item{metaFitSamRerun}{For batman rerun, a matrix (\eqn{t \times (s*n)}) containing 
 \eqn{s} posterior samples of total metabolites fit in its columns. \eqn{n} is the number 
 of spectra analyzed and \eqn{t} is the number of data points in each spectrum. 
 The remaining \eqn{n-1} spectra metabolites fit results are saved in the same sequence 
 in subsequent columns.}
 
 \item{metaIndFitSamRerun}{For batman rerun, a matrix (\eqn{t \times (m*s*n)}) 
 containing \eqn{s} posterior samples of \eqn{m} individual metabolites fit in its columns.
 \eqn{n} is the number of spectra analyzed and \eqn{t} is the number of data points 
 in each spectrum. The remaining \eqn{n-1} spectra results are saved 
 in the same sequence in subsequent columns.}
 
 \item{thetaSamRerun}{For batman rerun, a matrix (\eqn{t \times (s*n)}) containing 
 \eqn{s} samples of wavelet fit in its columns. \eqn{n} is the number of spectra analyzed.
 The remaining \eqn{n-1} spectra wavelet fit results are saved in the same sequence 
 in subsequent columns.}
 
 \item{outputDir}{The directory of output folder with all the output result files.}
}
\seealso{
\code{\link{batman}, \link{readBatmanOutput}}
}
\examples{
library(batman)
## Run batman
if(interactive()) 
{
bm<-batman()
## then call batmanrerun
bm<-batmanrerun(bm)
}
#######################################################################
## The following is an example of what will be displayed in R 
## and what value the user could input:
#######################################################################
## Rerunning batman for  500  iterations.
## percentage completed...                                                                          
##   |                                                                  |   0%
## Size of each spectrum is 382.
## Size of metabolite list is 22.
## Constructing chain data structure...
## time used is 1 seconds.
## Running MCMC...                                                                      
##   |==================================================================| 100%
## time used is 65 seconds.
## saving posteriors...
## 
## For rerun, time elapsed 
##  65.96 seconds.
## Reading in saved data in folder
## .../user_specified_dir/runBATMAN/BatmanOutput/07_Dec_17_35_53
## Completed.
########################################################################
## Alternatively if more than 1 spectrum are included without using fixed  
## effect, user will be asked to input whether to parallelize the analysis 
## between spectra.
########################################################################
## How many parallel processes (multicores) do you want to run the 
## multi-spectra analysis?
## (Enter 1 for running them sequentially.)
##
## Parallel processing of multi spectra currently cannot display 
## the progress bar (or any words), if you input is > 1, please be patient 
## for the results :)
##
## 1: 2  ## user input
##
## For rerun, time elapsed 
##   64.4 
## Reading in saved data in folder
## .../user_specified_dir/runBATMAN/BatmanOutput/07_Dec_17_35_53
## Completed.
#######################################################################
}
\keyword{datasets}

