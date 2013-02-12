// written by Dr. Jie Hao, Dr William Astle
#include "spec_template.h"


void spec_template::sample_log_fwhh()
{
    double log_fwhh_prop;
    double logratio = 0;  
    
    log_fwhh_prop=my_normrnd(log_fwhh_draw, sqrt(pars.log_fwhh_prop_var), rng);

    matrix Lprop (L);

    vector<double> decon(pars.p,0.0);

    for (int t = 0; t < pars.l; t++)
    { 
      FTems[t].build_curve(pars.x,log_fwhh_prop, pars.freq);

      if(WAVELETS)
	  {
    	  dwt(FTems[t].curve_prop,pars.nlev, pars.h_vec, decon); 
    	  Lprop[t] = decon;
	  }
      else
	  {
		  Lprop[t] = FTems[t].curve_prop;
	  }
    }
   
    logratio = calculate_metropolis_ratio_eta_local(&Lprop);
    
    double logprior=0.0;
    logprior+=normlogpdf(log_fwhh_prop, pars.log_fwhh_prior_mean, sqrt(pars.log_fwhh_prior_var));
    logprior-=normlogpdf(log_fwhh_draw, pars.log_fwhh_prior_mean, sqrt(pars.log_fwhh_prior_var));

    logratio+=logprior;
    
    //logproposal
    double logacceptance=log(my_rand(rng));
    if(logratio>logacceptance)
    {
        for (int t = 0; t < pars.l; t++)
		{
	       FTems[t].prop_accept_log_fwhh();
		}
        shape_accepcount=shape_accepcount+1;
        L=Lprop;
        log_fwhh_draw=log_fwhh_prop;
    }
}
    
