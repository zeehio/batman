\name{Batman-Output}
\alias{Batman-Output}
\title{BATMAN Outputs are Explained Here}
\description{\code{batman} and \code{batmanrerun} return the results as a data list with the objects described 
 in their individual function. They also put results in .txt format in a folder named 
 after the start execution time (\eqn{date\_month\_hours\_mins\_seconds}) within either folder ".../runBATMAN/BatmanOutput" 
 or folder "extdata" depending on \code{batman} input \code{createDir} settings.}
\value{
 \code{batman} and \code{batmanrerun} save their results in the following files in the output folder:  
  
 \item{beta_\eqn{i}_rr_\eqn{j}.txt}{A column vector (\eqn{m \times 1}) containing the estimated 
 posterior mean of relative concentrations for \eqn{m} fitted metabolites of spectrum \eqn{i}.  
 For \code{batman} results, \eqn{j} 
 is 0, and for \code{batmanrerun} results, \eqn{j} is 1.}
 
 \item{beta_sam_\eqn{i}_rr_\eqn{j}.txt}{A matrix (\eqn{m \times s}) with each row containing the \eqn{s} posterior samples of 
 the relative concentrations for one fitted metabolite of spectrum \eqn{i}. 
 \eqn{m} is the total number of fitted metabolites. For \code{batman} results, \eqn{j} 
 is 0, and for \code{batmanrerun} results, \eqn{j} is 1.} 
 
 \item{delta_draw_mean_\eqn{i}.txt}{A column vector (\eqn{M \times 1}) containing the posterior
 mean of \eqn{M} multiplets ppm shift from the pre-set ppm position value in 
 multi_data.csv or multi_data_user.csv of spectrum \eqn{i}.} 
 
 \item{delta_sam_\eqn{i}.txt}{A matrix (\eqn{s \times M}) containing the posterior samples 
 of \eqn{M} multiplets ppm shift. Every column correspond the shift posterior samples 
 of one multiplet for spectrum \eqn{i}. \eqn{M} is the sum of all 
 multiplets in the fitted metabolites.}
 
 \item{L_\eqn{i}.txt}{A matrix (\eqn{t \times M}) with each column as the template of one 
 fitted metabolite for spectrum \eqn{i} before fitting. \eqn{t}
 is the number of data points in each spectrum.}
  
 \item{lambda_sam_\eqn{i}_rr_\eqn{j}.txt}{A column vector (\eqn{s \times 1}) 
 containing \eqn{s} posterior samples of \eqn{\lambda} (a scalar global precision  
 parameter) for spectrum \eqn{i}. For \code{batman} results, \eqn{j} is 0, 
 and for \code{batmanrerun} results, \eqn{j} is 1.}
 
 \item{metabolitesListUsed.txt}{A column vector (\eqn{m \times 1}) containing the 
 \eqn{m} metabolite names 
 which have multiplets in/near the ppm region specified in batmanOptions.txt and used in the fitting.}
 
 \item{metaFit_sam_\eqn{i}_rr_\eqn{j}.txt}{A matrix (\eqn{t \times s}) containing 
 \eqn{s} posterior samples of total metabolites fit during MCMC iterations in 
 its columns for spectrum \eqn{i}. \eqn{t} is the number of data points in each spectrum. 
 For \code{batman} results, \eqn{j} is 0, and for \code{batmanrerun} results, \eqn{j} is 1.}
 
 \item{metaIndFit_sam_\eqn{i}_rr_\eqn{j}.txt}{A matrix (\eqn{t \times (m*s)}) 
 containing \eqn{s} posterior samples of \eqn{m} individual metabolites fit in 
 its columns for spectrum \eqn{i}. \eqn{t} is the number of data points in each 
 spectrum. Every \eqn{m} columns are the \eqn{m} individual metabolite fit 
 samples for one posterior sample. For \code{batman} results, \eqn{j} 
 is 0, and for \code{batmanrerun} results, \eqn{j} is 1.}
 
 \item{metaTemp_\eqn{i}_rr_\eqn{j}.txt}{A matrix \eqn{t \times m}) containing the posterior means of 
 \eqn{m} fitted metabolite templates in its columns (down sampled) after burn in
 for spectra \eqn{i}. \eqn{t} is the number of data points in each spectrum.
 For \code{batman} results, \eqn{j} is 0, and for \code{batmanrerun} 
 results, \eqn{j} is 1.}
 
 \item{metaTempHR_\eqn{i}_rr_\eqn{j}.txt}{A matrix (\eqn{t \times m}) containing the posterior means of 
 \eqn{m} fitted metabolite templates in its columns (without down sample) after burn in
 for spectra \eqn{i}. \eqn{t} is the number of data points (without down sample) 
 in each spectrum. For \code{batman} results, \eqn{j} 
 is 0, and for \code{batmanrerun} results, \eqn{j} is 1.}
 
 \item{MultipletsPpmShifts.txt}{A table (\eqn{M \times n}) containing the posterior means
 of multiplets ppm shift for \eqn{M} multiplets as its rows. \eqn{M} is the sum of all 
 multiplets in the fitted metabolites and \eqn{n} is the number of spectra analyzed.}
 
 \item{NMRdata_mod_\eqn{i}.txt}{A matrix (\eqn{t \times 2}) containing the input 
 spectrum \eqn{i} in its original resolution. The first column is ppm value, 
 and the second column is the \eqn{i}th spectrum intensity.}
 
 \item{RelCon.txt}{A table (\eqn{m \times n}) of the posterior means of relative 
 concentrations for \eqn{m} fitted metabolites and \eqn{n} spectra.} 
 
 \item{RelConCreInt.txt}{A table (\eqn{m \times 2n}) containing the 95\% credible 
 intervals (2.5\% and 97.5\%) for the relative concentrations of \eqn{m} 
 fitted metabolites for \eqn{n} spectra.} 
 
 \item{specFit_\eqn{i}_rr_\eqn{j}.txt}{A matrix (\eqn{t \times 5}) of BATMAN fit 
 results with five columns as:
 \deqn{[ppm, Original spectrum, Metabolites fit, Wavelet fit, Overall sum]} of spectrum \eqn{i}. 
 For \code{batman} results, \eqn{j} is 0, and for \code{batmanrerun} results, \eqn{j} is 1.}
  
 \item{specFitHR_\eqn{i}_rr_\eqn{j}.txt}{A column vector (\eqn{t \times 1}) of 
 metabolite fit result in the original resolution for spectrum \eqn{i}. \eqn{t} is the 
 number of data points (without down sample) in each spectrum.
 For \code{batman} results, \eqn{j} is 0,
 and for \code{batmanrerun} results, \eqn{j} is 1.}
 
 \item{theta_sam_\eqn{i}_rr_\eqn{j}.txt}{A matrix (\eqn{t \times s}) containing 
 \eqn{s} samples of wavelet fit during MCMC iterations in its columns for spectrum \eqn{i}. 
 For \code{batman} results, \eqn{j} is 0, and for \code{batmanrerun} results, \eqn{j} is 1.}
 
 \item{batmanOptions.txt}{The same file copied from batman input. This file will 
 be used by \code{batmanrerun}.}
 \item{metabolitesList.txt}{The same file copied from batman input.}
 \item{NMRdata.txt}{The same file copied from batman input.}
 If any plotting is performed, pdf files of the figure will be saved. For details, 
 please refer to each plotting functions.
}
