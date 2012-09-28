\name{Batman-Input}
\alias{Batman-Input}
\title{BATMAN Input Files are Explained Here}
\description{\code{batman} gets input parameters and metabolite templates information 
 from the input files explained here. The input files are in either 
 folder ".../runBATMAN/BatmanInput" or folder "extdata" depending on 
 \code{batman} arguments. 
 The user can modify the parameter values in the following input files (
 do not change the name of these files): 
        batmanOptions.txt, 
        metabolitesList.csv,
        multi_data.csv, 
        multi_data_user.csv,
        NMRdata.txt.}
\arguments{
 \item{batmanOptions.txt}{Option file to be used by \code{batman}. A copy of this 
 file in the output directory is used for \code{batmanrerun}.
 The parameters in \code{batmanOptions.txt} file are explained here with example
 input values. The parameters have to be listed in the particular order given here.
 
 

 
 
  \code{Include ppm ranges for analysis: (1.2, 1.6) (2.1, 2.8)}
  
  
   \itemize{
     \item{Put each set of ppm range in a pair of parentheses in the same line, 
      separate start and end ppm values with a comma, separate each set of 
      ppm range with space.}}
  
  
  \code{No. of spectra to be included: 1}
  
  
  \itemize{
     \item{ Integer, if no. > 1 and fixed effect (same concentration for all spectra) 
      is 0, user will be asked to choose 
      whether to parallelize fittings between spectra when running  
      \code{batman} or \code{rerunbatman}.}}
  
  
  \code{Lower limit for spectrum intensity: -0.5}
  
  
  \itemize{
     \item{Spectrum intensity smaller than the lower limit 
      will be replaced by the lower limit.}}
  
  
  \code{Normalisation factor: 20000}
  
  
\itemize{
     \item{The whole spectrum will be divided by the normalisation factor.}}
  
  
  \code{Down sampling: 3}
  
  
\itemize{
     \item{Integer, number of spectra variable will be reduced by the factor of 
      the input parameter, 3, in this case. For the example shown, the spectra variables
      with the index \eqn{1:3:end} will be used for analysis.}}


  \code{Save metabolites fit same as the original spectrum resolution (1/0): 1}


\itemize{
     \item{ Whether to save the metabolites fitting result in
      the original resolution without down sampling. 
      Input 1 for yes, and 0 for no.}}


  \code{Set seed for random number generation: 25}

\itemize{
     \item{ Random number generation seed, integer.}}
  
    
  \code{Stop burn in at iteration: 4000}


\itemize{
     \item{Integer, this is the number of burn in iterations. The number of iterations 
      after burn in will be asked when running \code{batman}. The posterior samples
      will be saved in the frequency specified by the next parameter.}}
  
  
  \code{Save results in every ? iterations: 5}
  
  \itemize{
     \item{Integer, save posterior samples for every 5 iterations.}}


  \code{Same concentration for all spectra (fixed effect) (1/0): 0}
  
  
\itemize{
     \item{Whether all the input spectra have the same metabolite concentrations (e.g.
      technical replicates). Input 1 for yes, and 0 for no.}}
      
      
  \code{Rerun iterations: 5000}
  
  
 \itemize{
     \item{Integer, this is the number of iterations for \code{batmanrerun}.
      The rerun will use fixed multiplets positions obtained 
      from running \code{batman}. There is no burnin for batman rerun.}}
      

  \code{Start temperature: 1000}
  
  
\itemize{
     \item{Sets the start temperature parameter of the likelihood of tempering. Higher temperature
      may need more burnin iterations to cool down.}}
      
          
  \code{Spectrometer frequency (MHz): 600}
  
  
\itemize{
     \item{Spectrometer used to collect the spectrum.}}
  

  \code{a: 0.00001}
      
      
  \code{b: 0.000000001}


\itemize{
     \item{Hyper parameters for the global precision priors 
      (\eqn{\lambda \sim Gamma(a,b/2)}) on wavelet coefficients.}}
       
        
  \code{Mean of the prior on mu: 0}
  
  
  \code{Variance of the prior on mu: 0.1}
  
  
  \code{Proposal variance for the Metropolis-Hastings sampler for mu: 0.002}
  
  
  \code{Variance of each of the priors on the nu_m: 0.0025}
    
    
  \code{Proposal variance for the Metropolis-Hastings sampler for each nu_m: 0.0001}
  
  
\itemize{
     \item{For peak width, \eqn{\gamma}, in ppm of metabolite \eqn{m}, the model for 
      \eqn{\gamma} is \eqn{\ln(\gamma)= \mu + \nu_m}
      where \eqn{\mu} is the spectrum wide average log-peakwidth and 
      \eqn{\nu_m} is a random effect on metabolite deviation from \eqn{\mu}.
      The mean of each prior on \eqn{\nu_m} is 0. Set the variance of the prior on 
      \eqn{\nu_m} to 0 to turn off the random effect on peak width to keep peaks at
      the same width. The user can keep the proposal variance parameters unchanged
      for most of the case.
    }}

  \code{mean of the prior on tau: -0.01}
  
\itemize{
     \item{Hyper priors (\eqn{\tau}) on negative wavelet coefficient (truncated normal).
      A more negative value means the wavelet fit will have more negative component.}}
  
 
  \code{steep: 2}
 
 
\itemize{
     \item{This parameter is inversely proportional to the variance of the prior on \eqn{\tau}.}}


  \code{rdelta: 0.030}
  
\itemize{
     \item{Prior of the truncation on ppm shift for all multiplets, individual prior 
  for each multiplet can be changed in the "multi_data.csv" file. Increase this 
  parameter to allow multiplets to shift more.}}}

 
 \item{metabolitesList.csv}{List of metabolite names to be fitted. Put "\%" in front 
 of the metabolite name to comment out any metabolite for batman analysis.}


 \item{multi_data.csv}{Multiplet template parameters file, obtained from the online 
 Human Metabolome Database (HMDB) version 2.5. The user can modify the parameters in 
 the template file and specify ppm positions, and normal distribution truncation of ppm
 shift parameters (a positive value applied as +/- on the distribution). 
 
 
 The columns are:
  
  
  \code{Metabolite}: The name of metabolite the multiplets belongs to. 	
  
  
  \code{pos_in_ppm}: The ppm position of the multiplets.
  
  
  \code{couple_code}: Coupling code.
  
  
  \code{J_constant}: J constant.
  
  
  \code{no_of_protons}: Number of protons in each multiplet.
  
  
  \code{overwrite_pos}: The default is "n" for not overwrite position, and in that case
  the value in "pos_in_ppm" is used for each multiplet. If user want to use a different 
  value from "pos_in_ppm", it should be put in this column.
                
  
  \code{overwrite_truncation}: The default is "n", and the default truncation value is obtained from 
  the user input truncation on ppm shift (rdelta) in \code{batmanOptions.txt}. 
  If the user wants to use different truncations for specific multiplets, it should be put in this column.
  This value will be used to calculate the ppm shift variance value (truncation/5)
  for the corresponding multiplets.
  
  
  \code{Include_multiplet}: The default is "1" and all multiplets belong to the listed 
  metabolites will be used. Set to "0" to exclude certain multiplet(s) from listed metabolite(s).}

 \item{multi_data_user.csv}{Metabolite template parameters file for user to add new metabolites in the
 same format as \code{multi_data.csv}. }

 \item{NMRdata.txt}{The file has ppm value as its first column, and real part of the 
 NMR spectrum in each of the subsequent columns. This file will be used when none of the input data argument is given.}
}
