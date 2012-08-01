#include "spec_template.h"


void spec_template::sample_log_fwhh_re(size_t t)
{
    double log_fwhh_re_prop;
    double logratio = 0;

    log_fwhh_re_prop=my_normrnd(FTems[t].log_fwhh_re_draw, sqrt(pars.log_fwhh_re_prop_var), rng);

    matrix Lprop (L);

    vector<double> decon(pars.p,0.0);

    FTems[t].build_curve_fwhh_re_prop(pars.x, log_fwhh_re_prop, log_fwhh_draw, pars.freq);

    if(WAVELETS)
	{
	    dwt(FTems[t].curve_prop,pars.nlev, pars.h_vec, decon); 
        Lprop[t] = decon;
	}
    else
    	Lprop[t] = FTems[t].curve_prop;        

    logratio = calculate_metropolis_ratio_eta_local(&Lprop);
    
    double logprior=0.0;
  
    logprior+=normlogpdf(log_fwhh_re_prop, 0.0, sqrt(pars.log_fwhh_re_prior_var));
    logprior-=normlogpdf(FTems[t].log_fwhh_re_draw, 0.0, sqrt(pars.log_fwhh_re_prior_var));
    
    //logratio
    //proposal
    //cancels out
    logratio+=logprior;
  
    //logproposal
    double logacceptance=log(my_rand(rng));

    if(logratio>logacceptance)
    {
      FTems[t].prop_accept_log_fwhh_re(log_fwhh_re_prop);
      shape_accepcount=shape_accepcount+1;
      L[t]=Lprop[t];
      FTems[t].log_fwhh_re_draw=log_fwhh_re_prop;
    }
    
}
    
