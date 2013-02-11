\name{batman-package}
\alias{batman-package}
\docType{package}
\title{
Bayesian AuTomated Metabolite Analyser for NMR spectra (BATMAN)
}

\description{
BATMAN deconvolves resonance peaks from NMR spectra of complex mixtures and obtains 
concentration estimates for the corresponding metabolites automatically.
This is achieved through a database of spectral profiles for known metabolites and 
a Bayesian Markov Chain Monte Carlo algorithm. Users have the options to  
specify the multiplet ppm position, position shift range, peak width range and so on.
Parallel processing is available if processing several spectra.
\tabular{ll}{\tab
\if{html}{\figure{batmanlogo.png}}
\if{latex}{\figure{batmanlogo.png}{options: width=1.5in}} \cr
}
}
\details{	
\tabular{ll}{
Package: \tab batman\cr
Type: \tab Package\cr
Version: \tab 1.1.0\cr
Date: \tab 2013-02-05\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}}
\author{
 Jie Hao <j.hao@ic.ac.uk>, William Astle <william.astle@mail.mcgill.ca>,   
 Maria De Iorio <m.deiorio@ucl.ac.uk>, Timothy Ebbels <t.ebbels@ic.ac.uk>
}
\references{
Astle W et al. (2012) A Bayesian Model of NMR Spectra for the Deconvolution and 
Quantification of Metabolites in Complex Biological Mixtures. Journal of the American Statistical Association,
Vol. 107, Iss. 500, 2012
http://arxiv.org/PS_cache/arxiv/pdf/1105/1105.2204v1.pdf


Hao J et al. (2012) BATMAN-an R package for the automated quantification of 
metabolites from NMR spectra using a Bayesian Model. BIOINFORMATICS, 2012:28(15):2088-2090.


http://arxiv.org/ftp/arxiv/papers/1112/1112.5794.pdf
}
\keyword{ package }

