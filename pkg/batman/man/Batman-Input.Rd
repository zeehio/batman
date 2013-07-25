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
        NMRdata.txt.
        }
\arguments{
 \item{batmanOptions.txt}{Option file to be used by \code{batman}. A copy of this 
 file in the output directory is used for \code{batmanrerun}.
 The parameters in \code{batmanOptions.txt} file are explained here with example
 input values. The parameters have to be listed in the particular order given here, 
 and do not leave empty lines in between except beginning with the comment character "\%".
 Please note that for version 1.0.9 and later, one more input line,
 
 "Use specified chemical shift for spectra (chemShiftperSpectra.csv) file (1/0): 0",
 
 is added at the end of this file.
 For earlier version users updated to this version, running \code{batman} will 
 add the above input line at the end of the file if missing.


 
 
  \code{ppmRange - ppm ranges for analysis: (1.2, 1.6) (2.1, 2.8)}
  
  
   \itemize{
     \item{Put each set of ppm range in a pair of parentheses in the same line, 
      separate start and end ppm values with a comma, separate each set of 
      ppm range with space. Note that, very small number of spectra variables
      may cause error in wavelet analysis, do not give very narrow ppm ranges and also check the "Down sampling:" 
      factor below, which used together, may also left very small number of spectra variables. }}
  
  
  \code{specNo - Ranges of spectra number to be included (e.g. 1,3-4 etc.): 1-3, 5}
  
  
  \itemize{
     \item{Integer, if no. > 1 and fixed effect (same concentration for all spectra) 
      is 0, user will be asked to choose 
      whether to parallelize fittings between spectra when running  
      \code{batman} or \code{rerunbatman}.}}
  
  
  \code{negThresh - Truncation threshold for negative intensities: -0.5}
  
  
  \itemize{
     \item{Spectrum intensity smaller than the lower limit 
      will be replaced by the lower limit.}}
  
  
  \code{scaleFac - Intensity scale factor: 20000}
  
  
\itemize{
     \item{The whole spectrum will be divided by the normalisation factor.}}
  
  
  \code{downSamp - Down sampling factor: 3}
  
  
\itemize{
     \item{Integer, number of spectra variable will be reduced by the factor of 
      the input parameter, 3, in this case. For the example shown, the spectra variables
      with the index \eqn{1:3:end} will be used for analysis.}}


  \code{hiresFlag - Save metabolite fit at resolution of original spectrum?}
  
  \code{(Yes - 1 / No - 0): 1}


\itemize{
     \item{Whether to save the metabolites fitting result in
      the original resolution without down sampling. 
      Input 1 for yes, and 0 for no.}}


  \code{randSeed  - Random number seed: 25}

\itemize{
     \item{Random number generation seed, integer.}}
  
    
  \code{nItBurnin - Number of burn-in iterations: 4000}


\itemize{
     \item{Integer, this is the number of burn in iterations. The number of iterations 
      after burn in will be asked when running \code{batman}. The posterior samples
      will be saved in the frequency specified by the next parameter. If changing 
      the range of spectrum causing fitting results inconsistent, this indicates that the burn in 
      stage hasn't found the best chemical shift. User may need to increase burn in iterations or reduce 
      prior truncation on ppm shift for each multiple (adjust parameter "rdelta" below). }}
  
  
  \code{thinning - Save MCMC state in every ? iterations: 5}
  
  \itemize{
     \item{Integer, save posterior samples for every 5 iterations.}}


  \code{cfeFlag - Same concentration for all spectra (fixed effect)?}
  
   \code{(Yes - 1 / No - 0): 0}
  
  
\itemize{
     \item{Whether all the input spectra have the same metabolite concentrations (e.g.
      technical replicates). Input 1 for yes, and 0 for no.}}
      
      
  \code{nItRerun - Number of iterations for batmanrerun: 5000}
  
  
 \itemize{
     \item{Integer, this is the number of iterations for \code{batmanrerun}.
      The rerun will use fixed multiplets positions obtained 
      from running \code{batman}. There is no burnin for batman rerun.}}
      

  \code{startTemp - Start temperature: 1000}
  
  
\itemize{
     \item{Sets the start temperature parameter of the likelihood of tempering. Higher temperature
      may need more burnin iterations to cool down.}}
      
          
  \code{specFreq - Spectrometer frequency (MHz): 600}
  
  
\itemize{
     \item{Spectrometer used to collect the spectrum.}}
  

  \code{a - Gamma-distributed with shape a: 0.00001}
      
      
  \code{b - Gamma-distributed with scale b: 0.000000001}


\itemize{
     \item{Hyper parameters for the global precision priors 
      (\eqn{\lambda \sim Gamma(a,b/2)}) on wavelet coefficients.}}
       
        
  \code{muMean - Mean of prior on global peak width (mu) in ln(Hz): 0}
  
  
  \code{muVar - Variance of prior on global peak width (mu) in ln(Hz): 0.1}
  
  
  \code{muVar_prop - Variance of proposal distribution for mu in ln(Hz): 0.002}
  
  
  \code{nuMVar - Variance of prior on peak width offset (nu_m) in ln(Hz): 0.0025}
    
    
  \code{nuMVarProp - Variance of proposal distribution for nu_m in ln(Hz): 0.0001}
  
  
\itemize{
     \item{For peak width, \eqn{\gamma}, in Hz of metabolite \eqn{m}, the model for 
      \eqn{\gamma} is \eqn{\ln(\gamma)= \mu + \nu_m}
      where \eqn{\mu} is the spectrum wide average log-peakwidth and 
      \eqn{\nu_m} is a random effect on metabolite deviation from \eqn{\mu}.
      The mean of each prior on \eqn{\nu_m} is 0. Set the variance of the prior on 
      \eqn{\nu_m} to 0 to turn off the random effect on peak width to keep peaks at
      the same width. The user can keep the proposal variance parameters unchanged
      for most of the case.
    }}

  \code{tauMean - mean of the prior on tau: -0.01}
  
\itemize{
     \item{Hyper priors (\eqn{\tau}) on negative wavelet coefficient (truncated normal).
      A more negative value means the wavelet fit will have more negative component.}}
  
 
  \code{tauPrec - inverse of variance of prior on tau: 2}
 
 
\itemize{
     \item{This parameter is inversely proportional to the variance of the prior on \eqn{\tau}.}}


  \code{rdelta - Truncation of the prior on peak shift (ppm): 0.030}
  
\itemize{
     \item{Prior of the truncation on ppm shift for all multiplets, individual prior 
  for each multiplet can be changed in the "multi_data.csv" file. Increase this 
  parameter to allow multiplets to shift more. Please note, increasing this value may need more burn in
  iterations to find the best chemical shift for multiplets.}}
  
   \code{csFlag - Specify chemical shift for each multiplet in each spectrum?}
   
    \code{(chemShiftperSpectra.csv file) (Yes - 1 / No - 0): 0}
  
\itemize{
     \item{Input "1" to use file "chemShiftperSpectra.csv" to specify chemical shift per multiplet and per spectrum. 
   Input "0" will not use that file. User can use the MATLAB tool "SplineFitBATMAN" provided to get more accurate 
   chemical shift per spectra for each multiplet. This tool will save  chemical shift information into "chemShiftperSpectra.csv".}}
  }

 
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
  0 = singlet, 
  1 = doublet, 
  2 = triplet,
3 = quartet,
4 = quintet,
5 = sextet,
6 = septet,
1,1 = doublet of doublets, 
1,2 = doublet of triplets,
2,1 = triplet of doublets,
2,2 = triplet of triplets,
1,3 = doublet of quartets,
3,1 = quartet of doublets, 
1,1,1 = doublet of doublet of doublets,
2,3 = triplet of quartets, 
3,2 = quartet of triplets, 
3,3 = quartet of quartets. 
If "-1" is inputted here, a user specified multiplet can be created. 
	An example can be found in file "multi_data_user.csv". If "-2" is inputted here, a multiplet with range
  specified in ppm in the field "J_constant" is used. Examples can be found in file "multi_data_user.csv".
 
  
  \code{J_constant}: J constant. If "-1" is inputted in the previous field "couple_code", J_constant/f are the offset
 	of peaks from the mutiplet position (f is the magnet frequency). Note that the spectra are shown in reverse ppm axis, so 
   a positive offset means peak at higher ppm value, and a negative offset is peak at lower ppm value. If "-2" is inputted 
   in the previous field "couple_code", the field here requires a two values input seperated by comma, which specifies the ppm
   range of the multiplet in the pure spectrum. Note in this case, the field "Metabolite" name will also be the .txt file name
   containing the pure spectrum (refer to \code{\link{createPureSpectraTemplate}}).
  
  
  \code{no_of_protons}: Number of protons in each multiplet. If "-1" is inputted for "couple_code", one or more 
	(corresponding to J_constant) values can be given here as peak weights. The sum of no_of_protons is 
	the number of protons in this multiplet. 
  
  
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
