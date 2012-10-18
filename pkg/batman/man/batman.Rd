\name{batman}
\alias{batman}
\title{Perform BATMAN and Plot Analysis Result}
\description{The main function, it performs metabolite and wavelet fitting to input NMR spectra, 
 plots fitting results, posterior distributions for relateive concentrations and peak positions, 
 and saves output. If the input \code{createDir = TRUE}, a folder name "runBATMAN" 
 will be created in specified directory, within which, two folders "BatmanInput" and 
 "BatmanOutput" are created. "BatmanInput" contains the input data files copied from 
 installed package folder "extdata". The user only needs to modify files in this
 folder to change the settings for running \code{batman}. The \code{batman} 
 output files are saved in "BatmanOutput" subfolders.}
\usage{
 batman(BrukerDataDir, txtFile, rData, createDir = TRUE, 
        runBATMANDir = getwd(), overwriteDir = FALSE, 
        figBatmanFit = TRUE, listMeta = FALSE, 
        figRelCon = FALSE, figMetaFit = FALSE)}
\arguments{
\item{BrukerDataDir}{The directory of the folder containing 1D Bruker spectral data files. 
 If not specified, spectral data will be read in from one of the following inputs prioritized 
 in the order: \code{txtFile}, \code{rData} and NMRdata.txt in "BatmanInput" folder.}
\item{txtFile}{The .txt file containing spectral data in the format of first column ppm, and 
 the second column the real part of spectrum. If not specified, spectral data will be read in 
 from one of the following inputs prioritized in the order: \code{BrukerDataDir}, 
 \code{rData} and NMRdata.txt in "BatmanInput" folder.}
\item{rData}{The R data file containing spectral data in the format of first column ppm, 
 and the second column the real part of spectrum. If not specified, spectral data will be read in 
 from one of the following inputs prioritized in the order: \code{BrukerDataDir}, 
 \code{txtFile} and NMRdata.txt in "BatmanInput" folder.}
\item{createDir}{If set TRUE, a new BATMAN work directory will be created specified by \code{runBATMANDir}. 
 If set FALSE, batman input will be obtained from the "extdata" folder in batman package installation directory, 
 and the batman output files will also be put within this folder.  The default is TRUE.}
\item{runBATMANDir}{User specified BATMAN work directory, the default is current work directory. 
 It will only work when \code{createDir} is set TRUE.}
\item{overwriteDir}{If folder "runBATMAN" exists, set TRUE to overwrite folder. The default is FALSE.}
\item{figBatmanFit}{Plot metabolites and wavelets fit if set TRUE. The default is TRUE.}
\item{listMeta}{Individual metabolite fit will also be shown in the plot if set TRUE. The default is FALSE.}
\item{figRelCon}{Plot posterior samples of the relative concentration for fitted 
 metabolites with 95\% credible interval if set TRUE. The default is FALSE.}
\item{figMetaFit}{If set TRUE, plot the posterior mean of the metabolites fit with 95\% 
 credible interval. The default is FALSE.}
}
\value{	
 It returns a data list with the following objects:  
  
 \item{specTitle}{A matrix (\eqn{2 \times n}) containing the spectrum number in its
 first row and the corresponding title of the spectrum in its second row.}
 
 \item{sFit}{A matrix \eqn{t \times 5n} of BATMAN fit results (down sampled). For 1 spectrum, it 
 is a matrix with 5 columns: \deqn{[ppm, original spectrum, 
 metabolites fit, wavelets fit, overall fit].} The "overall fit" is the posterior mean of 
 the BATMAN fit results after MCMC burn in iterations. Certain numbers of burn in iterations are used
 at the beginning of an MCMC run for finding a good starting point.
 \eqn{n} is the number of spectra, and \eqn{t} is the number of data points in each spectrum.} 
 
 \item{sFitHR}{A matrix \eqn{t \times 3n} of BATMAN fit results in the original resolution (without 
 down sample). For 1 spectrum, it is a matrix with 3 columns:   
 \deqn{[ppm, original spectrum, metabolites fit].} \eqn{n} is the number of spectra, and
 \eqn{t} is the number of data points (without down sample) in each spectrum.} 
 
 \item{beta}{A matrix (\eqn{m \times n}) containing the posterior means of relative 
 concentrations for \eqn{m} fitted metabolites and \eqn{n} spectra after burn in.}
 
 \item{betaSam}{A matrix (\eqn{m \times (s*n)}) containing (for the first spectrum)
 \eqn{s} posterior samples of the relative concentrations in its rows.
 \eqn{m} is the number of fitted metabolites. \eqn{n} is the number of spectra analyzed. 
 The subsequent columns contain the same format of data for the rest \eqn{n-1} spectra.}
 
 \item{betaCI}{A matrix (\eqn{m \times 2n}) containing the 95\% credible interval of 
 the relative concentrations for \eqn{m} fitted metabolites. Every pair of columns is 
 for one spectrum.}
 
 \item{metaTemp}{A matrix (\eqn{t \times (m*n)}) containing the posterior means of 
 \eqn{m} fitted metabolite templates in its columns (down sampled) after burn in.
 \eqn{n} is the number of spectra analyzed and \eqn{t} is the number of data 
 points in each spectrum.}
 
 \item{metaTempHR}{A matrix (\eqn{t \times (m*n)}) containing the posterior means of 
 \eqn{m} fitted metabolite templates in its columns (without down sample) after burn in.
 \eqn{n} is the number of spectra analyzed and \eqn{t} is the number of data 
 points (without down sample) in each spectrum.}
     
 \item{metaFitSam}{A matrix (\eqn{t \times (s*n)}) containing \eqn{s} posterior samples 
 of total metabolites fit during MCMC iterations in its columns. \eqn{n} is the number 
 of spectra analyzed and \eqn{t} is the number of data points in each spectrum. 
 The remaining \eqn{n-1} spectra metabolites fit results are saved in the same sequence 
 in subsequent columns.}
 
 \item{metaIndFitSam}{A matrix (\eqn{t \times (m*s*n)}) containing \eqn{s} posterior 
 samples of \eqn{m} individual metabolites fit during MCMC iterations in its columns.
 \eqn{n} is the number of spectra analyzed and \eqn{t} is the number of data points 
 in each spectrum. The remaining \eqn{n-1} spectra results are saved 
 in the same sequence in subsequent columns.}
 
 \item{thetaSam}{A matrix (\eqn{t \times (s*n)}) containing \eqn{s} samples of wavelet
 fit during MCMC iterations in its columns. \eqn{n} is the number of spectra analyzed.
 The remaining \eqn{n-1} spectra wavelet fit results are saved in the same sequence 
 in subsequent columns.}
 
 \item{delta}{A matrix (\eqn{M \times n}) containing posterior means of \eqn{M} 
 multiplets ppm shift of fitted metabolites in its rows. \eqn{M} is the sum of all 
 multiplets in the fitted metabolts. Each column of the matrix 
 corresponds to one spectrum. If only 1 spectrum is analyzed, delta is a column vector.}
 
 \item{deltaSam}{A matrix (\eqn{s \times (M*n)}) containing the posterior samples 
 of multiplets ppm shift. Every \eqn{M} columns correspond the shift posterior samples 
 of \eqn{M} multiplets for one spectrum. \eqn{M} is the sum of all 
 multiplets in the fitted metabolts and \eqn{n} is the number of spectra analyzed.}
  
 \item{outputDir}{The directory of output folder with all the output result files.}
}
\seealso{
\code{\link{readBatmanOutput}, \link{batmanrerun}}
}
\examples{
library(batman)
## Run BATMAN
bm<-batman()
100
1
## This will create the folder "runBATMAN" in current working directory, 
## within the folder "runBATMAN", a subfolder "BatmanInput" contains all the
## input files batman uses. Users can modify "metabolitesList.csv", 
## "batmanOptions.txt" and so on to change the settings of batman.
## Please check "BatmanInput" for details on how to adjust input parameters.
########################################################################
## The following is an example of what will be displayed in R 
## and what value the user could input: 
########################################################################
## batman...
## Enter number of post-burn-in iterations (burn-in currently set to  
## 400  iterations):
## 1: 100   ## user input
##
## Enter a number of choice from the menu below:
##  
## 1: Include the default template of multiplets in multi_data.csv file only.
## 2: Include the user input template of multiplets in multi_data_user.csv file only.
## 3: Include both the above files.
##
## Selection: 1  ## user input
## Loading multi_data.csv...
## Percentage completed...
##   |                                                                  |   0%
## Size of each spectrum is 393.
## Size of metabolite list is 22.
## Constructing chain data structure...
## time used is 0 seconds.
## Running MCMC...                                                                   
##   |===================================================               |  80%
## time used for burnin is 76 seconds.                                                                         
##   |==================================================================| 100%
## time used is 95 seconds.
## saving posteriors...
##  
## time elapsed 
##   95.61 
##  second.
## Reading in saved data in folder
## .../user_specified_dir/runBATMAN/BatmanOutput/07_Dec_17_19_18
## Completed.
########################################################################
## Alternatively if more than 1 spectrum are included without using fixed effect
## (in batmanOptions.txt file, set 
## "Same concentration for all spectra (fixed effect) (1/0): 0"),  
## user will be asked to input the following parameter: 
########################################################################
## How many parallel processes (multicores) do you want to run 
## the multi-spectra analysis?
## (Enter 1 for running them sequentially.)
##
##  Parallel processing of multi spectra currently cannot display 
##  progress bar (or any words), if you input is > 1, please be patient 
##  for the results :)
##
## 1: 2   ## user input
## time elapsed 
##   78.79 
##  second.
## Reading in saved data in folder
## .../user_specified_dir/runBATMAN/BatmanOutput/07_Dec_17_35_53
## Completed.
########################################################################
}
\keyword{datasets}
